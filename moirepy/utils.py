import numpy as np


def get_rotation_matrix(theta_rad: float) -> np.ndarray:
    """
    Computes a 2D rotation matrix for a given angle.

    Args:
        theta_rad (float): The rotation angle in radians.

    Returns:
        np.ndarray: A 2x2 rotation matrix that rotates a point
            counterclockwise by `theta_rad`.

    ```python
    >>> get_rotation_matrix(np.pi/2)
    array([[ 6.123234e-17, -1.000000e+00],
           [ 1.000000e+00,  6.123234e-17]])
    ```
    """
    return np.array(
        [
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad),  np.cos(theta_rad)]
        ]
    )

def are_coeffs_integers(v1, v2, v3, tol=1e-8):
    """
    Solves a v1 + b v2 = v3 for scalars a and b using Cramer's Rule,
    and checks if both are integers (within a tolerance).
    
    v1, v2, v3 are lists or tuples or np arrays of length 2
    """
    a1, a2 = v1
    b1, b2 = v2
    c1, c2 = v3

    det = a1 * b2 - a2 * b1
    if abs(det) < tol:
        return False  # no unique solution

    a = (c1 * b2 - c2 * b1) / det
    b = (a1 * c2 - a2 * c1) / det

    ret = abs(a - round(a)) < tol and abs(b - round(b)) < tol
    return ret