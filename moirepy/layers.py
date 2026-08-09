import numpy as np
from typing import List, Tuple, Dict, Optional, Union

class Layer:
    def __init__(
        self,
        lv1: np.ndarray,
        lv2: np.ndarray,
        basis_points: List[Union[Tuple[float, float, str], Tuple[float, float, float, str]]],
        neighbours: Dict[str, np.ndarray],
        pbc: bool = False,
        orbitals: int = 1,
    ):
        """
        Base Layer class.
        basis_points should be a list of tuples: (x, y, z, type_label) or (x, y, type_label)
        """
        self.lv1 = np.array([lv1[0], lv1[1], 0.0], dtype=np.float64) if len(lv1) == 2 else np.array(lv1, dtype=np.float64)
        self.lv2 = np.array([lv2[0], lv2[1], 0.0], dtype=np.float64) if len(lv2) == 2 else np.array(lv2, dtype=np.float64)
        
        # Enforce coordinate convention
        if not np.isclose(self.lv1[1], 0.0, atol=1e-8):
            raise ValueError("lv1 must be along the x-axis (lv1[1] == 0)")
        if self.lv2[1] <= 0.0:
            raise ValueError("lv2 must have a positive y-component")
            
        self.pbc = pbc
        self.orbitals = orbitals
        
        # Parse 3D basis points
        self.basis_points = []
        self.basis_labels = []
        for bp in basis_points:
            # Handle 2D or 3D input gracefully
            if len(bp) == 3: # (x, y, label) -> assume z=0
                self.basis_points.append([float(bp[0]), float(bp[1]), 0.0])
                self.basis_labels.append(bp[2])
            else: # (x, y, z, label)
                self.basis_points.append([float(bp[0]), float(bp[1]), float(bp[2])])
                self.basis_labels.append(bp[3])
                
        self.basis_points = np.array(self.basis_points, dtype=np.float64)
        self.basis_labels = np.array(self.basis_labels, dtype=str)
        self._neighbours = neighbours
        
        # Rotations/Translations will be applied to points
        self._rot_m = np.eye(3, dtype=np.float64)
        self._translation = np.zeros(3, dtype=np.float64)
        
        # State
        self._points = None
        self._labels = None
        self._bigger_points = None
        self._bigger_labels = None
        self._bigger_to_primary = None
        self._owner_id = None

    def apply_strain(self, strain_tensor: np.ndarray) -> None:
        """
        Applies a 2x2 strain tensor to the primitive lattice vectors lv1 and lv2.
        """
        if self._owner_id is not None:
            from .utils import LatticeAlreadyFinalisedError
            raise LatticeAlreadyFinalisedError("Cannot apply strain after layer is added to a MoireStack.")
            
        strain_tensor = np.array(strain_tensor, dtype=np.float64)
        if strain_tensor.shape != (2, 2):
            raise ValueError("strain_tensor must be a 2x2 matrix")
            
        self.lv1[:2] = (np.eye(2) + strain_tensor) @ self.lv1[:2]
        self.lv2[:2] = (np.eye(2) + strain_tensor) @ self.lv2[:2]
        
        # Re-enforce convention if strain broke it
        angle = np.arctan2(self.lv1[1], self.lv1[0])
        if abs(angle) > 1e-8:
            from .utils import get_rotation_matrix
            rot = get_rotation_matrix(-angle)
            self.lv1 = rot @ self.lv1
            self.lv2 = rot @ self.lv2
            
            for i in range(len(self.basis_points)):
                self.basis_points[i] = rot @ self.basis_points[i]

    def perform_rotation_translation(self, rot_deg: float, translation: Tuple[float, float] = (0,0)):
        from .utils import get_rotation_matrix
        self._rot_m = get_rotation_matrix(np.radians(rot_deg))
        
        if len(translation) == 2:
            self._translation = np.array([translation[0], translation[1], 0.0], dtype=np.float64)
        else:
            self._translation = np.array(translation, dtype=np.float64)

    def generate_points(self, mlv1: np.ndarray, mlv2: np.ndarray, n1: int = 1, n2: int = 1):
        """
        Generates points in the moire supercell defined by mlv1 and mlv2,
        tiled by n1 and n2.
        """
        if self._owner_id is None:
            # We assign a dummy owner if called standalone (mostly for testing)
            self._owner_id = 0
            
        mlv1 = np.array(mlv1, dtype=np.float64)[:2] * n1
        mlv2 = np.array(mlv2, dtype=np.float64)[:2] * n2
        
        S = np.column_stack((mlv1, mlv2))
        Minv = np.linalg.inv(S)
        
        P = np.column_stack((self.lv1[:2], self.lv2[:2]))
        Pinv = np.linalg.inv(P)
        
        corners = np.array([
            [0.0, 0.0],
            mlv1,
            mlv2,
            mlv1 + mlv2
        ])
        
        # We need to find the inverse mapping considering rotation
        # So we test corners against the inverse rotation
        rot_2d = self._rot_m[:2, :2]
        inv_rot = np.linalg.inv(rot_2d)
        
        corners_unrot = (inv_rot @ (corners - self._translation[:2]).T).T
        
        # Project into primitive fractional coords
        corners_prim = (Pinv @ corners_unrot.T).T
        
        min_i = int(np.floor(np.min(corners_prim[:, 0]))) - 1
        max_i = int(np.ceil(np.max(corners_prim[:, 0]))) + 1
        min_j = int(np.floor(np.min(corners_prim[:, 1]))) - 1
        max_j = int(np.ceil(np.max(corners_prim[:, 1]))) + 1
        
        points = []
        labels = []
        seen = set()
        
        for i in range(min_i, max_i + 1):
            for j in range(min_j, max_j + 1):
                for k, bp in enumerate(self.basis_points):
                    raw = i * self.lv1 + j * self.lv2 + bp
                    pt = self._rot_m @ raw + self._translation
                    
                    # Fractional in supercell
                    frac = Minv @ pt[:2]
                    
                    # Canonical dedup
                    frac_c = np.round(frac % 1.0, decimals=6) % 1.0
                    
                    # Also need a unique key for the atom position
                    # To prevent identical atoms (from PBC boundaries mapping to 0) being added multiple times
                    key = (int(frac_c[0] * 1e6), int(frac_c[1] * 1e6), k)
                    
                    # We check if it falls strictly inside the 0..1 bounding box
                    # Since frac_c is positive, any point that Maps strictly inside [0, 1) belongs to the cell
                    if 0.0 <= frac[0] < 1.0 + 1e-6 and 0.0 <= frac[1] < 1.0 + 1e-6:
                        # Wait, due to edge floating point issues, standard checking is brittle.
                        # We just trust canonical dedup and check if the original frac was roughly inside
                        # To be safe and capture all boundary points perfectly:
                        if -1e-5 <= frac[0] < 1.0 + 1e-5 and -1e-5 <= frac[1] < 1.0 + 1e-5:
                            if key not in seen:
                                seen.add(key)
                                points.append((frac_c[0], frac_c[1], pt))
                                labels.append(self.basis_labels[k])
                                
        # Sort for deterministic ordering
        points_sorted = sorted(zip(points, labels), key=lambda x: (x[0][0], x[0][1], x[0][2][-1]))
        
        self._points = np.array([x[0][2] for x in points_sorted], dtype=np.float64)
        self._labels = np.array([x[1] for x in points_sorted], dtype=str)
        
        # Build bigger_points (3x3 tiling)
        bigger_pts = []
        bigger_lbls = []
        self._bigger_to_primary = []
        
        # Generate 9 tiles
        tile_idx = 0
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                shift = di * np.array([mlv1[0], mlv1[1], 0.0]) + dj * np.array([mlv2[0], mlv2[1], 0.0])
                for p_idx, pt in enumerate(self._points):
                    bigger_pts.append(pt + shift)
                    bigger_lbls.append(self._labels[p_idx])
                    self._bigger_to_primary.append(p_idx)
                    
        self._bigger_points = np.array(bigger_pts, dtype=np.float64)
        self._bigger_labels = np.array(bigger_lbls, dtype=str)

    def first_nearest_neighbours(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Builds a KDTree to find intra-layer nearest neighbours based on self._neighbours.
        Returns:
            primary_indices, bigger_indices, labels_i, labels_j
        """
        if self._points is None:
            raise RuntimeError("Must call generate_points first")
            
        from scipy.spatial import cKDTree
        
        tree = cKDTree(self._bigger_points)
        
        i_idx = []
        j_idx = []
        l_i = []
        l_j = []
        
        for p_idx, pt in enumerate(self._points):
            label = self._labels[p_idx]
            if label not in self._neighbours:
                continue
                
            deltas = self._neighbours[label]
            if len(deltas) == 0:
                continue
                
            rotated_deltas = (self._rot_m @ deltas.T).T
            targets = pt + rotated_deltas
            
            dists, indices = tree.query(targets, k=1, distance_upper_bound=1e-5)
            
            for d, b_idx in zip(dists, indices):
                if d < 1e-5:
                    i_idx.append(p_idx)
                    j_idx.append(b_idx)
                    l_i.append(label)
                    l_j.append(self._bigger_labels[b_idx])
                    
        return np.array(i_idx), np.array(j_idx), np.array(l_i), np.array(l_j)
