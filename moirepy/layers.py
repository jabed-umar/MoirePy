import numpy as np
from typing import Tuple  # , List
import matplotlib.pyplot as plt
import math
from scipy.spatial import KDTree
import json

from .utils import get_rotation_matrix, LatticeAlreadyFinalisedError
from . import moirepy_rust as rbck  # import rust backend


class Layer:  # parent class
    def __init__(
        self,
        lv1: np.ndarray,
        lv2: np.ndarray,
        basis_points: list,
        neighbours: dict,
        pbc: bool = False,
        study_proximity: int = 1
    ) -> None:
        
        if lv1[1] != 0 or lv2[1] < 0:
            raise ValueError(
                """lv1 was expected to be along the x-axis,
                and lv2 should have a +ve y component
                Please refer to the documentation for more information: https://jabed-umar.github.io/MoirePy/find_theta/
                """
            )
        
        self._rust_lattice = rbck.Layer(
            lv1=lv1.tolist(),
            lv2=lv2.tolist(),
            basis_data=[tuple(p) for p in basis_points],
            neighbours=neighbours,
            pbc=pbc,
            study_proximity=study_proximity
        )

    def perform_rotation_translation(self, rot: float, translation = None) -> None:
        if self._rust_lattice.points is not None:
            raise RuntimeError("Cannot perform rotation and translation after points have been generated.")

        if translation is None:
            translation = (0.0, 0.0)

        try:
            translation_tuple = tuple(map(float, translation))
            if len(translation_tuple) != 2:
                raise ValueError
        except (TypeError, ValueError):
            raise ValueError("Translation must be a sequence of two numbers (x, y).")

        self._rust_lattice.perform_rotation_translation(rot, translation_tuple)

    def generate_points(
            self,
            mlv1: np.ndarray,
            mlv2: np.ndarray,
            mln1: int = 1,
            mln2: int = 1,
    ) -> None:
        assert isinstance(mln1, int) and mln1 > 0, "mln1 must be a positive integer."
        assert isinstance(mln2, int) and mln2 > 0, "mln2 must be a positive integer."
        
        # Ensure inputs are numpy arrays with float type
        mlv1 = np.asarray(mlv1, dtype=float)
        mlv2 = np.asarray(mlv2, dtype=float)
        
        self._rust_lattice.generate_points(mlv1, mlv2, mln1, mln2)


    def first_nearest_neighbours(self, points=None, types=None):
        if points is None:
            points = self.points
        if types is None:
            types = self.point_types

        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError(f"points must be (N, 2), got {points.shape}")
        
        # If types are passed as names, Rust handles the string-to-ID mapping
        return self._rust_lattice.first_nearest_neighbours(points, types)

    def get_neighbors_within_radius(self, query_points, radius: float):
        query_points = np.asarray(query_points, dtype=float)
        if query_points.ndim != 2 or query_points.shape[1] != 2:
            raise ValueError(f"query_points must be (N, 2), got {query_points.shape}")
        
        return self._rust_lattice.get_neighbors_within_radius(query_points, float(radius))

    def plot_lattice(self, plot_connections: bool = True, colours: list = ["r", "g", "b", "c", "m", "y", "k"]) -> None:
        # Create a color map: index -> color
        # Since self.point_type_names gives ["A", "B", ...], we map 0 -> 'r', 1 -> 'g', etc.
        if len(colours) == 1: cols = {t[-1]:colours[0] for i, t, in enumerate(self.basis_points)}
        else: cols = {t[-1]:colours[i] for i, t, in enumerate(self.basis_points)}

        if plot_connections:
            a = []
            for this, point_type in zip(self.points, self.point_types):
                for delta in self.neighbours[point_type]:
                    a.append([this[0], this[0] + delta[0]])
                    a.append((this[1], this[1] + delta[1]))
                    a.append("k--")
            plt.plot(*a, alpha=0.1)

        # plt.title("Lattice Points")
        # plt.xlabel("X Coordinate")
        # plt.ylabel("Y Coordinate")
        plt.axis("equal")


    def __repr__(self):
        return (
            f"""Layer(
    lv1 = {self.lv1},
    lv2 = {self.lv2},
    basis_points = {self.basis_points},
    study_proximity = {self.study_proximity},
    pbc = {self.pbc},
)"""
        )


    # --- Read-Only Redirects to Rust ---
    @property
    def lv1(self):
        return self._rust_lattice.lv1
    @lv1.setter
    def lv1(self, value):
        raise LatticeAlreadyFinalisedError("lv1", self.__class__.__name__)

    @property
    def lv2(self):
        return self._rust_lattice.lv2
    @lv2.setter
    def lv2(self, value):
        raise LatticeAlreadyFinalisedError("lv2", self.__class__.__name__)

    @property
    def basis_points(self):
        return self._rust_lattice.basis_points
    @basis_points.setter
    def basis_points(self, value):
        raise LatticeAlreadyFinalisedError("basis_points", self.__class__.__name__)

    @property
    def basis_types(self):
        return self._rust_lattice.basis_types
    @basis_types.setter
    def basis_types(self, value):
        raise LatticeAlreadyFinalisedError("basis_types", self.__class__.__name__)

    @property
    def neighbours(self):
        return dict(zip(self.basis_types, self._rust_lattice.neighbours))
    @neighbours.setter
    def neighbours(self, value):
        raise LatticeAlreadyFinalisedError("neighbours", self.__class__.__name__)

    @property
    def mlv1(self):
        return self._rust_lattice.mlv1
    @mlv1.setter
    def mlv1(self, value):
        raise LatticeAlreadyFinalisedError("mlv1", self.__class__.__name__)

    @property
    def mlv2(self):
        return self._rust_lattice.mlv2
    @mlv2.setter
    def mlv2(self, value):
        raise LatticeAlreadyFinalisedError("mlv2", self.__class__.__name__)

    @property
    def mln1(self):
        return self._rust_lattice.mln1
    @mln1.setter
    def mln1(self, value):
        raise LatticeAlreadyFinalisedError("mln1", self.__class__.__name__)

    @property
    def mln2(self):
        return self._rust_lattice.mln2
    @mln2.setter
    def mln2(self, value):
        raise LatticeAlreadyFinalisedError("mln2", self.__class__.__name__)

    @property
    def pbc(self):
        return self._rust_lattice.pbc
    @pbc.setter
    def pbc(self, value):
        raise LatticeAlreadyFinalisedError("pbc", self.__class__.__name__)

    @property
    def study_proximity(self):
        return self._rust_lattice.study_proximity
    @study_proximity.setter
    def study_proximity(self, value):
        raise LatticeAlreadyFinalisedError("study_proximity", self.__class__.__name__)

    @property
    def points(self):
        return self._rust_lattice.points
    @points.setter
    def points(self, value):
        raise LatticeAlreadyFinalisedError("points", self.__class__.__name__)

    @property
    def point_types(self):
        return self._rust_lattice.point_types
    @point_types.setter
    def point_types(self, value):
        raise LatticeAlreadyFinalisedError("point_types", self.__class__.__name__)

    @property
    def rot_m(self):
        return self._rust_lattice.rot_m
    @rot_m.setter
    def rot_m(self, value):
        raise LatticeAlreadyFinalisedError("rot_m", self.__class__.__name__)

    @property
    def translation(self):
        return self._rust_lattice.translation
    @translation.setter
    def translation(self, value):
        raise LatticeAlreadyFinalisedError("translation", self.__class__.__name__)


# ===============================================
# ============= some example layers =============
# ===============================================

class SquareLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        # local definitions, passed to parent
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0, 1])  # Lattice vector in the y-direction
        basis_points = [
            # location of the point inside the unit cell
            (0, 0, "A"),
        ]
        neighbours = {
            "A": [
                [1, 0],   # Right
                [0, 1],   # Up
                [-1, 0],  # Left
                [0, -1],  # Down
            ],
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, study_proximity)



class TriangularLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3) / 2])  # Lattice vector at 60 degrees
        basis_points = [
            (0, 0, "A"),
        ]
        neighbours = {
            "A": [
                [1, 0],  # Right
                [0.5, np.sqrt(3) / 2],  # Right-up
                [-0.5, np.sqrt(3) / 2],  # Left-up
                [-1, 0],  # Left
                [-0.5, -np.sqrt(3) / 2],  # Left-down
                [0.5, -np.sqrt(3) / 2],  # Right-down
            ],
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, study_proximity)



class Rhombus60Layer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        angle = 60  # hardcoded angle... make a copy of the whole class for different angles
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        cos_angle = np.cos(np.radians(angle))
        sin_angle = np.sin(np.radians(angle))
        lv2 = np.array([cos_angle, sin_angle])  # Lattice vector at specified angle
        basis_points = [
            (0, 0, "A"),
        ]
        neighbours = {
            "A": [
                [1, 0],  # Right
                [cos_angle, sin_angle],  # Up
                [-1, 0],  # Left
                [-cos_angle, -sin_angle],  # Down
            ],
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, study_proximity)



class KagomeLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3)/2])  # Lattice vector at 60 degrees

        basis_points = [
            (0, 0, "A"),
            (0.5, 0, "B"),
            (0.25, np.sqrt(3)/4, "C"),
        ]

        neighbours = {
            "A": [
                [ 0.5,              0],  # Right
                [ 0.25,  np.sqrt(3)/4],  # Right-up
                [-0.5,              0],  # Left
                [-0.25, -np.sqrt(3)/4],  # Left-down
            ],
            "B": [
                [ 0.5,              0],  # Right
                [-0.25,  np.sqrt(3)/4],  # Left-up
                [-0.5,              0],  # Left
                [ 0.25, -np.sqrt(3)/4],  # Right-down
            ],
            "C": [
                [ 0.25,  np.sqrt(3)/4],  # Right-up
                [-0.25,  np.sqrt(3)/4],  # Left-up
                [-0.25, -np.sqrt(3)/4],  # Left-down
                [ 0.25, -np.sqrt(3)/4],  # Right-down
            ],
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, study_proximity)


class HexagonalLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0]) # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3) / 2])

        basis_points = [
            # coo_x, coo_y, point_type (unique)
            (0, 0, "A"),
            (1, 1/np.sqrt(3), "B"),
        ]
        neighbours = {
            "A": [
                [0, 1/np.sqrt(3)],
                [-0.5, -1/(2 * np.sqrt(3))],
                [ 0.5, -1/(2 * np.sqrt(3))],
            ],
            "B": [
                [0.5, 1/(2 * np.sqrt(3))],
                [-0.5, 1/(2 * np.sqrt(3))],
                [0, -1/np.sqrt(3)],
            ],
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, study_proximity)


