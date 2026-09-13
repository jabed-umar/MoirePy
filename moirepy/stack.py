import numpy as np
from typing import List, Optional, Tuple
from .utils import StackStateError
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree
import scipy.sparse as sp
from .layers import Layer

from . import moirepy_rust

class MoireStack:
    def __init__(self, base_layer: Layer):
        """
        Initializes a MoireStack with a base layer (Layer 0).
        Twists of subsequent layers are relative to this layer.
        """
        if base_layer._owner_id is not None:
            raise ValueError("Layer is already part of a MoireStack!")
            
        self.layers: List[Layer] = [base_layer]
        base_layer._owner_id = id(self)
        
        self.lv1_ref = base_layer.lv1[:2]
        self.lv2_ref = base_layer.lv2[:2]
        
        self.mlv1: Optional[np.ndarray] = None
        self.mlv2: Optional[np.ndarray] = None
        
        self.state = "BUILDING"

    def add_layer(
        self, 
        layer: Layer, 
        ll1: int, 
        ll2: int, 
        ul1: int, 
        ul2: int, 
        translate: Tuple[float, float] = (0.0, 0.0),
        z_pos: Optional[float] = None
    ):
        if len(translate) != 2:
            if len(translate) == 3:
                raise ValueError(
                    "The 'translate' argument only accepts a 2D vector (x, y) for in-plane shifts. "
                    "For out-of-plane (z) shifts, please use the explicit 'z_pos' argument instead. "
                    f"Example: translate=({translate[0]}, {translate[1]}), z_pos={translate[2]}"
                )
            else:
                raise ValueError("translate must be a 2D tuple (x, y).")
                
        import warnings
        if z_pos is None or np.isclose(z_pos, 0.0, atol=1e-8):
            warnings.warn(
                "z_pos was not provided or is exactly 0.0. "
                "This means the added layer will physically collide in the same 2D plane as the base layer. "
                "Did you forget to specify the interlayer distance (e.g., z_pos=3.35)?",
                UserWarning, stacklevel=2
            )
            if z_pos is None:
                z_pos = 0.0
        if self.state != "BUILDING":
            raise StackStateError("Cannot add layer after stack is finalized.")
            
        if layer._owner_id is not None:
            raise ValueError("Layer is already part of a MoireStack!")
            
        # First added layer defines the supercell
        if self.mlv1 is None:
            self.mlv1 = ll1 * self.lv1_ref + ll2 * self.lv2_ref
            
            # mlv2 is mlv1 rotated by the natural angle of the base lattice
            # Angle between lv1 and lv2
            cos_beta = np.dot(self.lv1_ref, self.lv2_ref) / (np.linalg.norm(self.lv1_ref) * np.linalg.norm(self.lv2_ref))
            beta = np.arccos(np.clip(cos_beta, -1.0, 1.0))
            
            # Since lv2[1] > 0 and lv1[1] == 0, beta is positive
            cos_b = np.cos(beta)
            sin_b = np.sin(beta)
            rot_beta = np.array([
                [cos_b, -sin_b],
                [sin_b, cos_b]
            ])
            self.mlv2 = rot_beta @ self.mlv1
            
        # Vector mapped by the new layer's indices
        new_vec = ul1 * layer.lv1[:2] + ul2 * layer.lv2[:2]
        
        # Strict Math Commensurability check!
        l_mlv1 = np.linalg.norm(self.mlv1)
        l_new = np.linalg.norm(new_vec)
        
        if not np.isclose(l_mlv1, l_new, atol=1e-6):
            raise ValueError(
                f"Commensurate length mismatch! mlv1 length: {l_mlv1:.6f}, "
                f"new layer vector length: {l_new:.6f}. "
                "Did you copy the correct values from the website? Or forget to apply_strain?"
            )
            
        # Calculate twist angle to align new_vec to mlv1
        # angle from new_vec to mlv1
        angle_mlv1 = np.arctan2(self.mlv1[1], self.mlv1[0])
        angle_new = np.arctan2(new_vec[1], new_vec[0])
        
        theta_rad = angle_mlv1 - angle_new
        theta_deg = np.degrees(theta_rad)
        
        layer.perform_rotation_translation(theta_deg, translate, z_pos)
        layer._owner_id = id(self)
        self.layers.append(layer)


    def generate_points(self, n1: int = 1, n2: int = 1):
        """
        Locks the stack geometry and generates points for all layers.
        """
        if self.state != "BUILDING":
            raise StackStateError("generate_points can only be called in BUILDING state.")
            
        if self.mlv1 is None:
            # Only one layer in stack?
            self.mlv1 = self.lv1_ref
            self.mlv2 = self.lv2_ref
            
        for layer in self.layers:
            layer.generate_points(self.mlv1, self.mlv2, n1, n2)
            
        self.state = "GENERATED"

    def repeat(self, n1: int, n2: int):
        """
        Tiles the supercell n1 x n2 times.
        Can only be called in GENERATED state.
        """
        if self.state != "GENERATED":
            raise StackStateError("repeat can only be called in GENERATED state.")
            
        for layer in self.layers:
            # Re-run layer generation with the scaled boundaries
            layer.generate_points(self.mlv1, self.mlv2, n1, n2)
            
        # Update the stack's supercell vectors to match the repeated size
        self.mlv1 = self.mlv1 * n1
        self.mlv2 = self.mlv2 * n2

    def generate_connections(self, radii_matrix: np.ndarray) -> None:
        """
        Builds interlayer connections up to given radii matrix.
        radii_matrix should be an L x L matrix where L is number of layers.
        """
        if self.state != "GENERATED":
            raise StackStateError("generate_connections must be called in GENERATED state.")
            
        radii_matrix = np.array(radii_matrix, dtype=np.float64)
        num_layers = len(self.layers)
        if radii_matrix.shape != (num_layers, num_layers):
            raise ValueError(f"radii_matrix must be {num_layers}x{num_layers}")
            
        
        self.connections = {} # (i, j) -> (idx_i, idx_j, label_i, label_j, R_shift)
        
        for i in range(num_layers):
            for j in range(num_layers):
                if i == j:
                    continue
                r = radii_matrix[i, j]
                if r <= 0:
                    continue
                    
                pts_i = np.ascontiguousarray(self.layers[i]._points, dtype=np.float64)
                bigger_pts_j = np.ascontiguousarray(self.layers[j]._bigger_points, dtype=np.float64)
                prim_pts_j = np.ascontiguousarray(self.layers[j]._points, dtype=np.float64)
                b2p_j = np.ascontiguousarray(self.layers[j]._bigger_to_primary, dtype=np.int64)
                
                c_i, c_j, R_shifts = moirepy_rust.find_neighbors_within_radius(
                    pts_i,
                    bigger_pts_j,
                    prim_pts_j,
                    b2p_j,
                    float(r)
                )
                
                l_i = self.layers[i]._labels[c_i]
                l_j = self.layers[j]._labels[c_j]
                
                self.connections[(i, j)] = (
                    c_i,
                    c_j,
                    l_i,
                    l_j,
                    R_shifts
                )

        
        trees = []
        for l in self.layers:
            trees.append(cKDTree(l._bigger_points))
            
        for i in range(num_layers):
            for j in range(num_layers):
                if i == j:
                    continue
                r = radii_matrix[i, j]
                if r <= 0:
                    continue
                    
                pts_i = self.layers[i]._points
                tree_j = trees[j]
                
                indices_list = tree_j.query_ball_point(pts_i, r)
                
                c_i = []
                c_j = []
                l_i = []
                l_j = []
                R_shifts = []
                
                for idx_i, n_indices in enumerate(indices_list):
                    for b_idx_j in n_indices:
                        idx_j = self.layers[j]._bigger_to_primary[b_idx_j]
                        
                        pt_i = self.layers[i]._points[idx_i]
                        b_pt_j = self.layers[j]._bigger_points[b_idx_j]
                        pt_j = self.layers[j]._points[idx_j]
                        
                        shift = b_pt_j - pt_j
                        
                        c_i.append(idx_i)
                        c_j.append(idx_j)
                        l_i.append(self.layers[i]._labels[idx_i])
                        l_j.append(self.layers[j]._bigger_labels[b_idx_j])
                        R_shifts.append(shift)
                        
                R_shifts_arr = np.array(R_shifts, dtype=np.float64)
                if len(R_shifts) == 0:
                    R_shifts_arr = np.zeros((0, 3), dtype=np.float64)
                    
                self.connections[(i, j)] = (
                    np.array(c_i, dtype=np.int64),
                    np.array(c_j, dtype=np.int64),
                    np.array(l_i, dtype=str),
                    np.array(l_j, dtype=str),
                    R_shifts_arr
                )

    def generate_hamiltonian(self, t_matrix: dict, t_self: dict) -> None:
        """
        Builds the Hamiltonian as a sparse COO matrix block by block.
        t_matrix: (i, j) -> callable(coo_i, coo_j, R, labels_i, labels_j, stack) (R is supercell shift)
        t_self: i -> callable(coo, labels, stack)
        """
        
        num_layers = len(self.layers)
        
        # Calculate block offsets based on orbitals
        offsets = [0]
        for l in self.layers:
            offsets.append(offsets[-1] + len(l._points) * l.orbitals)
        total_dim = offsets[-1]
            
        blocks_to_assemble = []
        
        for i in range(num_layers):
            for j in range(num_layers):
                if i == j:
                    if i in t_self:
                        func = t_self[i]
                        i_idx, j_idx, l_i, l_j = self.layers[i].first_nearest_neighbours()
                        if hasattr(self.layers[i], '_cached_intra_shifts') and self.layers[i]._cached_intra_shifts is not None:
                            j_primary, R_shifts = self.layers[i]._cached_intra_shifts
                        else:
                            b2p = self.layers[i]._bigger_to_primary
                            j_primary = b2p[j_idx]
                            R_shifts = np.ascontiguousarray(
                                self.layers[i]._bigger_points[j_idx] - self.layers[i]._points[j_primary],
                                dtype=np.float64
                            )
                            self.layers[i]._cached_intra_shifts = (j_primary, R_shifts)
                        
                        call_args = (i_idx, j_primary, R_shifts, l_i, l_j, self)
                        blocks_to_assemble.append((
                            call_args,
                            func,
                            i_idx,
                            j_primary,
                            R_shifts,
                            offsets[i],
                            offsets[i],
                            False
                        ))
                else:
                    pair = (i, j)
                    if pair in t_matrix or (j, i) in t_matrix:
                        if pair not in t_matrix:
                            # Auto-conjugate from (j, i)
                            func = t_matrix[(j, i)]
                            if (j, i) in self.connections:
                                c_j, c_i, l_j, l_i, R_shifts = self.connections[(j, i)]
                                call_args = (c_j, c_i, R_shifts, l_j, l_i, self)
                                blocks_to_assemble.append((
                                    call_args,
                                    func,
                                    c_i,
                                    c_j,
                                    R_shifts,
                                    offsets[i],
                                    offsets[j],
                                    True
                                ))
                        else:
                            func = t_matrix[pair]
                            if pair in self.connections:
                                c_i, c_j, l_i, l_j, R_shifts = self.connections[pair]
                                call_args = (c_i, c_j, R_shifts, l_i, l_j, self)
                                blocks_to_assemble.append((
                                    call_args,
                                    func,
                                    c_i,
                                    c_j,
                                    R_shifts,
                                    offsets[i],
                                    offsets[j],
                                    False
                                ))
                                
        h_data, h_row, h_col, rx_data, ry_data, rz_data = moirepy_rust.assemble_hamiltonian_sparse(blocks_to_assemble)
        
        self.H = coo_matrix((h_data, (h_row, h_col)), shape=(total_dim, total_dim))
        self.Rx = coo_matrix((rx_data, (h_row, h_col)), shape=(total_dim, total_dim))
        self.Ry = coo_matrix((ry_data, (h_row, h_col)), shape=(total_dim, total_dim))
        self.Rz = coo_matrix((rz_data, (h_row, h_col)), shape=(total_dim, total_dim))
            

    def get_phase(self, k: np.ndarray) -> np.ndarray:
        """
        Calculates the bloch phase factor e^{i k . (r_j - r_i + R_shift)} for each non-zero element in H.
        Since H only stores non-zero elements, we just need to return an array of same length as H.data.
        Ensure phase expands to match orbitals block dimensions for perfect H.multiply(phase) usage.
        """
        if not hasattr(self, 'H'):
            raise RuntimeError("Must call generate_hamiltonian first.")
            
        k = np.array(k, dtype=np.float64)
        if k.shape == (2,):
            k = np.array([k[0], k[1], 0.0])
            
        # Get absolute positions of all points
        pts = []
        for l in self.layers:
            for pt in l._points:
                for _ in range(l.orbitals):
                    pts.append(pt)
        all_pts = np.array(pts)
        
        # d = r_j - r_i + R_shift
        r_i = all_pts[self.H.row]
        r_j = all_pts[self.H.col]
        
        shift_x = self.Rx.data
        shift_y = self.Ry.data
        shift_z = self.Rz.data
        shifts = np.column_stack((shift_x, shift_y, shift_z))
        
        d = r_j - r_i + shifts
        
        # Phase = e^{i k . d}
        k_dot_d = np.sum(d * k, axis=1)
        phase = np.exp(1j * k_dot_d)
        
        return phase

    def plot_3d(self, show_connections: bool = True, enable_hover: bool = True):
        """
        Renders an interactable 3D Plotly graph of the MoireStack.
        Allows toggling individual layers and connection lines natively via the legend.
        Set enable_hover=True in Jupyter to dynamically highlight bonds.
        """
        try:
            import plotly.graph_objects as go
        except ImportError:
            raise ImportError("Plotly is required for 3D plotting. Install it via 'pip install plotly'.")
            
        if self.state == "BUILDING":
            raise StackStateError("Must call generate_points() before plotting.")
            
        if enable_hover:
            fig = go.FigureWidget()
        else:
            fig = go.Figure()
        # fig = go.FigureWidget()
        
        # Plot layers
        for i, layer in enumerate(self.layers):
            pts = layer._points
            if pts is None or len(pts) == 0:
                continue
                
            fig.add_trace(go.Scatter3d(
                x=pts[:, 0].tolist(),
                y=pts[:, 1].tolist(),
                z=pts[:, 2].tolist(),
                mode='markers',
                marker=dict(size=4, opacity=0.9),
                name=f"Layer {i}"
            ))
            
            # Draw internal bonds (honeycomb) natively for the layer
            i_idx, b_idx, l_i, l_j = layer.first_nearest_neighbours()
            if len(i_idx) > 0:
                pts_i = pts[i_idx]
                pts_j = layer._bigger_points[b_idx]
                
                x_lines = np.empty(len(pts_i) * 3, dtype=object)
                x_lines[0::3] = pts_i[:, 0]
                x_lines[1::3] = pts_j[:, 0]
                x_lines[2::3] = None
                
                y_lines = np.empty(len(pts_i) * 3, dtype=object)
                y_lines[0::3] = pts_i[:, 1]
                y_lines[1::3] = pts_j[:, 1]
                y_lines[2::3] = None
                
                z_lines = np.empty(len(pts_i) * 3, dtype=object)
                z_lines[0::3] = pts_i[:, 2]
                z_lines[1::3] = pts_j[:, 2]
                z_lines[2::3] = None
                
                fig.add_trace(go.Scatter3d(
                    x=x_lines.tolist(),
                    y=y_lines.tolist(),
                    z=z_lines.tolist(),
                    mode='lines',
                    line=dict(width=2),
                    opacity=0.5,
                    name=f"Bonds L{i}"
                ))
            
        # Plot connections
        if show_connections and hasattr(self, 'connections'):
            for (i, j), (c_i, c_j, l_i, l_j, R_shifts) in self.connections.items():
                if j <= i:
                    continue # Skip intralayer (i==j) and avoid duplicate reverse interlayer (j<i)
                    
                pts_i = self.layers[i]._points[c_i]
                pts_j = self.layers[j]._points[c_j] + R_shifts
                
                # To draw lines in Plotly rapidly, we interleave [x1, x2, None, x3, x4, None...]
                x_lines = np.empty(len(pts_i) * 3, dtype=object)
                x_lines[0::3] = pts_i[:, 0]
                x_lines[1::3] = pts_j[:, 0]
                x_lines[2::3] = None
                
                y_lines = np.empty(len(pts_i) * 3, dtype=object)
                y_lines[0::3] = pts_i[:, 1]
                y_lines[1::3] = pts_j[:, 1]
                y_lines[2::3] = None
                
                z_lines = np.empty(len(pts_i) * 3, dtype=object)
                z_lines[0::3] = pts_i[:, 2]
                z_lines[1::3] = pts_j[:, 2]
                z_lines[2::3] = None
                
                fig.add_trace(go.Scatter3d(
                    x=x_lines.tolist(),
                    y=y_lines.tolist(),
                    z=z_lines.tolist(),
                    mode='lines',
                    line=dict(color='gray', width=1),
                    opacity=0.3,
                    visible='legendonly',
                    name=f"Interlayer (L{i}-L{j})"
                ))
                
        # Plot lattice boundaries (Supercell) at each layer's z-height
        if self.mlv1 is not None and self.mlv2 is not None:
            v1 = self.mlv1
            v2 = self.mlv2
            bx = [0, v1[0], v1[0] + v2[0], v2[0], 0]
            by = [0, v1[1], v1[1] + v2[1], v2[1], 0]
            
            for i, layer in enumerate(self.layers):
                if layer._points is not None and len(layer._points) > 0:
                    z_val = np.mean(layer._points[:, 2])
                    # z_val = layer._translation[2]
                    bz = [z_val] * 5
                    
                    fig.add_trace(go.Scatter3d(
                        x=bx,
                        y=by,
                        z=bz,
                        mode='lines',
                        line=dict(color='black', width=3, dash='dash'),
                        name=f"Boundary L{i}",
                        showlegend=False
                    ))
                
        fig.update_layout(
            title="MoirePy 3D Lattice",
            scene=dict(
                xaxis_title='X (Å)',
                yaxis_title='Y (Å)',
                zaxis_title='Z (Å)',
                aspectmode='data', # Critical: forces physical aspect ratio!
                dragmode='orbit',  # Orbital rotation instead of turntable
                camera=dict(
                    projection=dict(type='orthographic'), # No perspective shrinking
                    eye=dict(x=0, y=0, z=2.5),            # Look from the top (along Z axis)
                    up=dict(x=0, y=1, z=0)                # Up is Y axis
                )
            ),
            margin=dict(l=0, r=0, b=0, t=40)
        )
        
        # --- INTERACTIVE HOVER HIGHLIGHTING ---
        if enable_hover:
            # Add a dedicated trace for fast dynamic highlighting in Jupyter
            fig.add_trace(go.Scatter3d(
                x=[], y=[], z=[],
                mode='lines',
                line=dict(color='red', width=5),
                name='Highlighted Bonds',
                hoverinfo='none'
            ))
    
            def highlight_bonds(trace, points, state):
                if not points.point_inds: return
                idx = points.point_inds[0]
                hovered_pt = np.array([trace.x[idx], trace.y[idx], trace.z[idx]], dtype=np.float64)
                
                hx, hy, hz = [], [], []
                for t in fig.data:
                    if t.mode == 'lines' and ('Interlayer' in t.name or 'Bonds' in t.name):
                        for i in range(0, len(t.x), 3):
                            if t.x[i] is None: 
                                continue
                            pt1 = np.array([t.x[i], t.y[i], t.z[i]], dtype=np.float64)
                            pt2 = np.array([t.x[i+1], t.y[i+1], t.z[i+1]], dtype=np.float64)
                            if np.allclose(pt1, hovered_pt) or np.allclose(pt2, hovered_pt):
                                hx.extend([pt1[0], pt2[0], None])
                                hy.extend([pt1[1], pt2[1], None])
                                hz.extend([pt1[2], pt2[2], None])
                                
                with fig.batch_update():
                    fig.data[-1].x = hx
                    fig.data[-1].y = hy
                    fig.data[-1].z = hz
                            
            for t in fig.data:
                if t.mode == 'markers':
                    t.on_hover(highlight_bonds)
        # ----------------------------------------

        return fig

