import numpy as np

def get_rotation_matrix(theta):
    theta_rad = np.deg2rad(theta)
    return np.array(
        [
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad),  np.cos(theta_rad)]
        ]
    )
