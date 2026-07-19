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
