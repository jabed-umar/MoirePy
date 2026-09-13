import numpy as np
from .layers import Layer

class SquareLayer(Layer):
    def __init__(self, pbc: bool = False, orbitals: int = 1):
        lv1 = np.array([1.0, 0.0, 0.0])
        lv2 = np.array([0.0, 1.0, 0.0])
        basis_points = [(0.0, 0.0, 0.0, "A")]
        neighbours = {
            "A": np.array([
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0]
            ])
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, orbitals)

class TriangularLayer(Layer):
    def __init__(self, pbc: bool = False, orbitals: int = 1):
        lv1 = np.array([1.0, 0.0, 0.0])
        lv2 = np.array([0.5, np.sqrt(3.0)/2.0, 0.0])
        basis_points = [(0.0, 0.0, 0.0, "A")]
        neighbours = {
            "A": np.array([
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.5, np.sqrt(3.0)/2.0, 0.0],
                [-0.5, -np.sqrt(3.0)/2.0, 0.0],
                [-0.5, np.sqrt(3.0)/2.0, 0.0],
                [0.5, -np.sqrt(3.0)/2.0, 0.0]
            ])
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, orbitals)

class HexagonalLayer(Layer):
    def __init__(self, pbc: bool = False, orbitals: int = 1):
        lv1 = np.array([1.0, 0.0, 0.0])
        lv2 = np.array([0.5, np.sqrt(3.0)/2.0, 0.0])
        basis_points = [
            (0.0, 0.0, 0.0, "A"),
            (0.5, 0.5/np.sqrt(3.0), 0.0, "B")
        ]
        neigh_deltas = np.array([
            [0.0, 1.0/np.sqrt(3.0), 0.0],
            [-0.5, -0.5/np.sqrt(3.0), 0.0],
            [0.5, -0.5/np.sqrt(3.0), 0.0]
        ])
        neighbours = {
            "A": -neigh_deltas,
            "B": neigh_deltas
        }
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, orbitals)

class KagomeLayer(Layer):
    def __init__(self, pbc: bool = False, orbitals: int = 1):
        lv1 = np.array([1.0, 0.0, 0.0])
        lv2 = np.array([0.5, np.sqrt(3.0)/2.0, 0.0])
        basis_points = [
            (0.0, 0.0, 0.0, "A"),
            (0.5, 0.0, 0.0, "B"),
            (0.25, np.sqrt(3.0)/4.0, 0.0, "C")
        ]
        # Neighbours simplified for testing
        neighbours = {"A": [], "B": [], "C": []}
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, orbitals)

class Rhombus60Layer(Layer):
    def __init__(self, pbc: bool = False, orbitals: int = 1):
        lv1 = np.array([1.0, 0.0, 0.0])
        lv2 = np.array([0.5, np.sqrt(3.0)/2.0, 0.0])
        basis_points = [(0.0, 0.0, 0.0, "A")]
        neighbours = {"A": []}
        super().__init__(lv1, lv2, basis_points, neighbours, pbc, orbitals)
