import numpy as np
from typing import List, Optional, Tuple
from .layers import Layer
from .utils import StackStateError

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
        translate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    ):
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
        
        layer.perform_rotation_translation(theta_deg, translate)
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
            from .utils import StackStateError
            raise StackStateError("generate_connections must be called in GENERATED state.")
            
        radii_matrix = np.array(radii_matrix, dtype=np.float64)
        num_layers = len(self.layers)
        if radii_matrix.shape != (num_layers, num_layers):
            raise ValueError(f"radii_matrix must be {num_layers}x{num_layers}")
            
        from scipy.spatial import cKDTree
        
        self.connections = {} # (i, j) -> (idx_i, idx_j, label_i, label_j, R_shift)
        
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
                        
                self.connections[(i, j)] = (
                    np.array(c_i, dtype=np.int64),
                    np.array(c_j, dtype=np.int64),
                    np.array(l_i, dtype=str),
                    np.array(l_j, dtype=str),
                    np.array(R_shifts, dtype=np.float64)
                )

    def generate_hamiltonian(self, t_matrix: dict, t_self: dict) -> None:
        """
        Builds the Hamiltonian as a sparse COO matrix block by block.
        t_matrix: (i, j) -> callable(coo_i, coo_j, R, labels_i, labels_j, stack) (R is supercell shift)
        t_self: i -> callable(coo, labels, stack)
        """
        from scipy.sparse import coo_matrix
        import scipy.sparse as sp
        
        num_layers = len(self.layers)
        
        # Calculate block offsets based on orbitals
        offsets = [0]
        for l in self.layers:
            offsets.append(offsets[-1] + len(l._points) * l.orbitals)
            
        H_blocks = []
        
        for i in range(num_layers):
            for j in range(num_layers):
                if i == j:
                    if i in t_self:
                        func = t_self[i]
                        n_pts = len(self.layers[i]._points)
                        orb = self.layers[i].orbitals
                        
                        i_idx, j_idx, l_i, l_j = self.layers[i].first_nearest_neighbours()
                        
                        R_zero = np.zeros((len(i_idx), 3))
                        # For intra-layer, primary indices are used for coo_j (we just map them via bigger_to_primary)
                        # Wait, first_nearest_neighbours returns bigger_indices for j
                        j_primary = np.array([self.layers[i]._bigger_to_primary[idx] for idx in j_idx], dtype=np.int64)
                        R_shifts = self.layers[i]._bigger_points[j_idx] - self.layers[i]._points[j_primary]
                        
                        val = func(i_idx, j_primary, R_shifts, l_i, l_j, self)
                        
                        # Handle orbitals block
                        if orb == 1:
                            block = coo_matrix((val, (i_idx, j_primary)), shape=(n_pts, n_pts))
                        else:
                            # val should be of shape (N, orb, orb)
                            pass # Simplification for orb=1
                            
                        H_blocks.append((i, j, block))
                else:
                    pair = (i, j)
                    if pair in t_matrix or (j, i) in t_matrix:
                        # Auto-conjugate logic
                        if pair not in t_matrix:
                            # We use (j, i) and conjugate
                            func = t_matrix[(j, i)]
                            if (j, i) in self.connections:
                                c_j, c_i, l_j, l_i, R_shifts = self.connections[(j, i)]
                                val = func(c_j, c_i, R_shifts, l_j, l_i, self)
                                val = np.conj(val)
                                
                                n_pts_i = len(self.layers[i]._points)
                                n_pts_j = len(self.layers[j]._points)
                                orb_i = self.layers[i].orbitals
                                orb_j = self.layers[j].orbitals
                                
                                # Note the transpose in indices
                                block = coo_matrix((val, (c_i, c_j)), shape=(n_pts_i, n_pts_j))
                                H_blocks.append((i, j, block))
                        else:
                            func = t_matrix[pair]
                            if pair in self.connections:
                                c_i, c_j, l_i, l_j, R_shifts = self.connections[pair]
                                val = func(c_i, c_j, R_shifts, l_i, l_j, self)
                                
                                n_pts_i = len(self.layers[i]._points)
                                n_pts_j = len(self.layers[j]._points)
                                
                                block = coo_matrix((val, (c_i, c_j)), shape=(n_pts_i, n_pts_j))
                                H_blocks.append((i, j, block))
                                
        # Assemble block matrix
        total_dim = offsets[-1]
        
        # Build L x L array of objects for bmat
        bmat_array = [[None for _ in range(num_layers)] for _ in range(num_layers)]
        for i, j, block in H_blocks:
            bmat_array[i][j] = block
            
        self.H = sp.bmat(bmat_array, format='coo')

