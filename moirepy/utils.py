import numpy as np

class LatticeAlreadyFinalisedError(Exception):
    """Raised when attempting to modify a lattice that is already finalised."""
    pass

class StackStateError(Exception):
    """Raised when an operation is performed in the wrong MoireStack phase."""
    pass

def get_rotation_matrix(theta_rad: float) -> np.ndarray:
    """
    Returns a 3x3 rotation matrix for rotating a 3D vector around the Z-axis by theta_rad.
    The Z coordinate remains invariant.
    """
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    return np.array([
        [cos_t, -sin_t, 0.0],
        [sin_t,  cos_t, 0.0],
        [0.0,    0.0,   1.0]
    ], dtype=np.float64)

def are_coeffs_integers(v1: np.ndarray, v2: np.ndarray, v3: np.ndarray, tol: float = 1e-8) -> bool:
    """
    Checks if v3 is an integer linear combination of v1 and v2.
    Expects 2D or 3D vectors. Only the first two components (x, y) are used for the check,
    since the periodicity is 2D.
    """
    # Use only the x and y components
    basis = np.column_stack((v1[:2], v2[:2]))
    try:
        coeffs = np.linalg.solve(basis, v3[:2])
    except np.linalg.LinAlgError:
        return False
    
    # Check if coefficients are integers within tolerance
    return bool(np.all(np.abs(coeffs - np.round(coeffs)) < tol))
