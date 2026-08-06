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
