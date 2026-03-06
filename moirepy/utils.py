import numpy as np


def get_rotation_matrix(theta_rad: float) -> np.ndarray:
    """Return the 2D counterclockwise rotation matrix.

    Parameters
    ----------
    theta_rad : float
        Rotation angle in radians.

    Returns
    -------
    np.ndarray
        Rotation matrix of shape ``(2, 2)``.

    Examples
    --------
    >>> get_rotation_matrix(np.pi / 2)
    array([[ 6.123234e-17, -1.000000e+00],
           [ 1.000000e+00,  6.123234e-17]])
    """
    return np.array(
        [
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad),  np.cos(theta_rad)]
        ]
    )

def are_coeffs_integers(v1, v2, v3, tol=1e-8):
    """Check if ``v3`` is an integer linear combination of ``v1`` and ``v2``.

    Solves ``a * v1 + b * v2 = v3`` in 2D using Cramer's rule, then tests
    whether ``a`` and ``b`` are integers up to a tolerance.

    Parameters
    ----------
    v1, v2, v3 : array-like of shape (2,)
        Input 2D vectors.
    tol : float, optional
        Tolerance for determinant singularity and integer closeness test.

    Returns
    -------
    bool
        ``True`` if a unique solution exists and both coefficients are
        integer-valued within ``tol``; otherwise ``False``.
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

# class LatticeAlreadyFinalisedError(RuntimeError):
class LatticeAlreadyFinalisedError(Exception):
    """Error raised when attempting to mutate immutable lattice parameters.

    Parameters
    ----------
    varname : str
        Name of the attempted-to-modify attribute.
    classname : str
        Class name of the immutable object instance.
    """

    def __init__(self, varname: str, classname: str):
        message = (
            f"Variable '{varname}' is not editable after creation."
            f"Please make a new instance of {classname} instead."
        )
        super().__init__(message)
