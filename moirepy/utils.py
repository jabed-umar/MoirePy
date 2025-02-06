import numpy as np


def get_rotation_matrix(theta_rad: float) -> np.array:
    """
    Computes a 2D rotation matrix for a given angle.

    Args:
        theta_rad (float): The rotation angle in radians.

    Returns:
        np.array: A 2x2 rotation matrix that rotates a point
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
