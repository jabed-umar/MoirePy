import numpy as np
from typing import List, Tuple
import matplotlib.pyplot as plt

from utils import get_rotation_matrix



class Layer:
    def __init__(self) -> None:
        self.toll = self.lattice_constant * 1e-3
        rot_m = get_rotation_matrix(self.rotation)
        
        # Rotate lv1 and lv2 vectors
        self.lv1 = rot_m @ self.lv1
        self.lv2 = rot_m @ self.lv2

        # Rotate lattice_points
        self.lattice_points = [
            [*(rot_m @ np.array([x, y])), atom_type]
            for x, y, atom_type in self.lattice_points
        ]

        # Rotate neighbours
        self.neighbours = {
            atom_type: [rot_m @ np.array(neighbour) for neighbour in neighbour_list]
            for atom_type, neighbour_list in self.neighbours.items()
        }

        self.lattice_points: List[Tuple[float, float]] = self.generate_points()

    def generate_points(self) -> dict:
        points = {atom_type: [] for _, _, atom_type in self.lattice_points}

        for i in range(self.nx):
            for j in range(self.ny):
                for x, y, atom_type in self.lattice_points:
                    new_point = i * self.lv1 + j * self.lv2 + np.array([x, y])
                    points[atom_type].append(new_point)
                    new_point = i * self.lv1 - j * self.lv2 + np.array([x, y])
                    points[atom_type].append(new_point)
                    new_point = -i * self.lv1 + j * self.lv2 + np.array([x, y])
                    points[atom_type].append(new_point)
                    new_point = -i * self.lv1 - j * self.lv2 + np.array([x, y])
                    points[atom_type].append(new_point)

        return points

    def plot_lattice(self, plot_connections: bool = True, plot_unit_cell: bool = False) -> None:
        # plt.figure(figsize=(8, 8))

        for atom_type, atom_points in self.lattice_points.items():
            x_coords = [point[0] for point in atom_points]
            y_coords = [point[1] for point in atom_points]
            plt.scatter(x_coords, y_coords, s=50)

            if plot_connections:
                for point in atom_points:
                    for neighbor in self.neighbours[atom_type]:
                        connection = point + np.array(neighbor)
                        plt.plot(
                            [point[0], connection[0]],
                            [point[1], connection[1]],
                            "r--",
                            alpha=0.5,
                        )

        if plot_unit_cell:
            for i in range(self.ny + 1):
                # line from (lv1*0 + lv2*i) to (lv1*nx + lv2*i)
                plt.plot(
                    [self.lv1[0] * 0 + self.lv2[0] * i, self.lv1[0] * self.nx + self.lv2[0] * i],
                    [self.lv1[1] * 0 + self.lv2[1] * i, self.lv1[1] * self.nx + self.lv2[1] * i],
                    "k:",
                    alpha=0.3,
                )
            
            for i in range(self.nx + 1):
                # line from (lv1*i + lv2*0) to (lv1*i + lv2*ny)
                plt.plot(
                    [self.lv1[0] * i + self.lv2[0] * 0, self.lv1[0] * i + self.lv2[0] * self.ny],
                    [self.lv1[1] * i + self.lv2[1] * 0, self.lv1[1] * i + self.lv2[1] * self.ny],
                    "k:",
                    alpha=0.3,
                )

        plt.title("Lattice Points")
        plt.xlabel("X Coordinate")
        plt.ylabel("Y Coordinate")
        plt.axis("equal")




class HexagonalLayer(Layer):
    def __init__(self, rotation: float = 0, lattice_constant=1) -> None:
        self.rotation = rotation
        self.lattice_constant = lattice_constant
        self.lv1 = np.array([np.sqrt(3), 0]) * lattice_constant
        self.lv2 = np.array([np.sqrt(3) / 2, 3 / 2]) * lattice_constant
        
        self.lattice_points = (
            # coo_x, coo_y, atom_type (unique)
            [0, 0, "A"],
            [np.sqrt(3), 1, "B"],
        )
        self.neighbours = {
            "A": [
                [0, 1],
                [-np.sqrt(3) / 2, -0.5],
                [np.sqrt(3) / 2, -0.5],
            ],
            "B": [
                [np.sqrt(3) / 2, 0.5],
                [-np.sqrt(3) / 2, 0.5],
                [0, -1],
            ],
        }
        super().__init__()


class SquareLayer(Layer):
    def __init__(self, nx: int, ny: int, rotation: float = 0, lattice_constant: float | int = 1) -> None:
        self.nx = nx
        self.ny = ny
        self.rotation = rotation
        self.lattice_constant = lattice_constant
        self.lv1 = np.array([1, 0]) * lattice_constant  # Lattice vector in the x-direction
        self.lv2 = np.array([0, 1]) * lattice_constant  # Lattice vector in the y-direction
        self.lattice_points = (
            [0, 0, "A"],
        )
        self.neighbours = {
            "A": [
                [1, 0],  # Right
                [0, 1],  # Up
                [-1, 0], # Left
                [0, -1], # Down
            ],
        }
        super().__init__()

class TriangularLayer(Layer):
    def __init__(self, nx: int, ny: int, rotation: float = 0, lattice_constant: float | int = 1) -> None:
        self.nx = nx
        self.ny = ny
        self.rotation = rotation
        self.lattice_constant = lattice_constant
        self.lv1 = np.array([1, 0]) * lattice_constant  # Lattice vector in the x-direction
        self.lv2 = np.array([0.5, np.sqrt(3)/2]) * lattice_constant  # Lattice vector at 60 degrees
        self.lattice_points = (
            [0, 0, "A"],
        )
        self.neighbours = {
            "A": [
                [1, 0],  # Right
                [0.5, np.sqrt(3)/2],  # Right-up
                [-0.5, np.sqrt(3)/2],  # Left-up
                [-1, 0],  # Left
                [-0.5, -np.sqrt(3)/2],  # Left-down
                [0.5, -np.sqrt(3)/2],  # Right-down
            ],
        }
        super().__init__()

class RhombusLayer(Layer):
    def __init__(self, nx: int, ny: int, rotation: float = 0, angle: float = 60, lattice_constant: float | int = 1, restrict: bool = True) -> None:
        self.nx = nx
        self.ny = ny
        self.rotation = rotation
        self.lattice_constant = lattice_constant
        
        if restrict and (120 < angle or angle < 60):
            raise ValueError("Angle must be between 60 and 120 degrees (inclusive). "
                             "If this is not true, the neighbours will not be calculated correctly. "
                             "However, you can set restrict=False to bypass this check.")
        
        self.lv1 = np.array([1, 0]) * lattice_constant  # Lattice vector in the x-direction
        cos_angle = np.cos(np.radians(angle))
        sin_angle = np.sin(np.radians(angle))
        self.lv2 = np.array([cos_angle, sin_angle]) * lattice_constant  # Lattice vector at specified angle
        self.lattice_points = (
            [0, 0, "A"],
        )
        self.neighbours = {
            "A": [
                [1, 0],  # Right
                [cos_angle, sin_angle],  # Up
                [-1, 0],  # Left
                [-cos_angle, -sin_angle],  # Down
            ],
        }
        super().__init__()


class KagomeLayer(Layer):
    def __init__(self, nx: int, ny: int, rotation: float = 0, lattice_constant: float | int = 1) -> None:
        self.nx = nx
        self.ny = ny
        self.rotation = rotation
        self.lattice_constant = lattice_constant
        self.lv1 = np.array([1, 0]) * 2 * lattice_constant  # Lattice vector in the x-direction
        self.lv2 = np.array([0.5, np.sqrt(3)/2]) * 2 * lattice_constant  # Lattice vector at 60 degrees

        self.lattice_points = (
            [0, 0, "A"],
            [1, 0, "B"],
            [0.5, np.sqrt(3)/2, "C"],
        )

        self.neighbours = {
            "A": [
                [1, 0],  # Right
                [0.5, np.sqrt(3)/2],  # Right-up
                [-1, 0],  # Left
                [-0.5, -np.sqrt(3)/2],  # Left-down
            ],
            "B": [
                [1, 0],  # Right
                [-0.5, np.sqrt(3)/2],  # Left-up
                [-1, 0],  # Left
                [0.5, -np.sqrt(3)/2],  # Right-down
            ],
            "C": [
                [0.5, np.sqrt(3)/2],  # Right-up
                [-0.5, np.sqrt(3)/2],  # Left-up
                [-0.5, -np.sqrt(3)/2],  # Left-down
                [0.5, -np.sqrt(3)/2],  # Right-down
            ],
        }

        super().__init__()



if __name__ == "__main__":
    # layer1 = HexagonalLayer (10, 10, 0)
    # layer1 = SquareLayer    (20, 20, 0)
    # layer1 = TriangularLayer(10, 10, 0)
    layer1 = RhombusLayer   (10, 10, 0)
    # layer1 = KagomeLayer    (10, 10, 0)
    
    # layer2 = HexagonalLayer (10, 10, 10)
    # layer2 = SquareLayer    (20, 20, 22.6198)
    # layer2 = TriangularLayer(10, 10, 10)
    # layer2 = RhombusLayer   (20, 20, 10)
    layer2 = KagomeLayer    (5, 5, 10)
    
    # print("initialized layer")
    
    layer1.plot_lattice(0, 0)
    layer2.plot_lattice(0, 0)
    plt.grid()
    plt.show()
