from typing import Union, Callable
import numpy as np
from layers import Layer, HexagonalLayer, SquareLayer, RhombusLayer, TriangularLayer, KagomeLayer
import matplotlib.pyplot as plt
from utils import get_rotation_matrix


class MoireLattice:
    def __init__(self, latticetype: Layer, a:int, b:int, nx:int=1, ny:int=1, pbc = False):
        lower_lattice = latticetype(pbc=pbc)
        upper_lattice = latticetype(pbc=pbc)

        lv1, lv2 = lower_lattice.lv1, lower_lattice.lv2

        # c = cos(theta) between lv1 and lv2
        c = np.dot(lv1, lv2) / (np.linalg.norm(lv1) * np.linalg.norm(lv2))
        beta = np.arccos(c)
        mlv1 = lv1 * a + lv2 * b
        mlv2 = get_rotation_matrix(beta).dot(mlv1)
        
        # the actual theta is the angle between a*lv1 + b*lv2 and b*lv1 + a*lv2
        one = a * lv1 + b * lv2
        two = b * lv1 + a * lv2
        c = np.dot(one, two) / (np.linalg.norm(one) * np.linalg.norm(two))
        theta = np.arccos(c)  # in radians
        print(f"theta = {theta:.4f} rad ({np.rad2deg(theta):.4f} deg)")

        upper_lattice.perform_rotation(theta)

        lower_lattice.generate_points(mlv1, mlv2, nx, ny)
        upper_lattice.generate_points(mlv1, mlv2, nx, ny)

        self.a = a
        self.b = b
        self.nx = nx
        self.ny = ny
        self.lower_lattice = lower_lattice
        self.upper_lattice = upper_lattice
        self.theta = theta
        self.mlv1 = mlv1
        self.mlv2 = mlv2
        self.ham = None

        print(f"{len(self.lower_lattice.points)} points in lower lattice")
        print(f"{len(self.upper_lattice.points)} points in upper lattice")

        # self.plot_lattice()

    def plot_lattice(self):
        mlv1 = self.mlv1
        mlv2 = self.mlv2
        nx = self.nx
        ny = self.ny

        plt.plot(*zip(*self.lower_lattice.points), 'r.')
        plt.plot(*zip(*self.upper_lattice.points), 'b.')

        # parallellogram around the whole lattice
        plt.plot([0, nx*mlv1[0]], [0, nx*mlv1[1]], 'k', linewidth=1)
        plt.plot([0, ny*mlv2[0]], [0, ny*mlv2[1]], 'k', linewidth=1)
        plt.plot([nx*mlv1[0], nx*mlv1[0] + ny*mlv2[0]], [nx*mlv1[1], nx*mlv1[1] + ny*mlv2[1]], 'k', linewidth=1)
        plt.plot([ny*mlv2[0], nx*mlv1[0] + ny*mlv2[0]], [ny*mlv2[1], nx*mlv1[1] + ny*mlv2[1]], 'k', linewidth=1)

        # just plot mlv1 and mlv2 parallellogram
        plt.plot([0, mlv1[0]], [0, mlv1[1]], 'k', linewidth=1)
        plt.plot([0, mlv2[0]], [0, mlv2[1]], 'k', linewidth=1)
        plt.plot([mlv1[0], mlv1[0] + mlv2[0]], [mlv1[1], mlv1[1] + mlv2[1]], 'k', linewidth=1)
        plt.plot([mlv2[0], mlv1[0] + mlv2[0]], [mlv2[1], mlv1[1] + mlv2[1]], 'k', linewidth=1)

        plt.grid()
        plt.show()

    def generate_hamiltonian(
        self,
        tuu: Union[float, int, Callable[[float], float]] = 1,
        tdd: Union[float, int, Callable[[float], float]] = 1,
        tdu: Union[float, int, Callable[[float], float]] = 1,
    ):
        if isinstance(tuu, int) or  isinstance(tuu, float): tuu = lambda theta: tuu
        if isinstance(tdd, int) or  isinstance(tdd, float): tdd = lambda theta: tdd
        if isinstance(tdu, int) or  isinstance(tdu, float): tdu = lambda theta: tdu
        assert callable(tuu) and callable(tdd) and callable(tdu), "tuu, tdd and tdu must be floats, ints or callable objects like functions"

        # hamiltonian has 4 parts:
        # 1. interation inside the lower lattice
        # 2. interaction inside the upper lattice
        # 3. interaction between the lower and upper lattice
        # 4. interaction between the upper and lower lattice
        
        # 1. interaction inside the lower lattice
        ham_ll = np.zeros((len(self.lower_lattice.points), len(self.lower_lattice.points)))
        # for 
        

# lattice = MoireLattice(TriangularLayer, 9, 10, 3+0, 2+0, pbc=True)
# lattice = MoireLattice(TriangularLayer, 2, 3, 1, 1, pbc=True)
lattice = MoireLattice(TriangularLayer, 9, 10, 2, 4, pbc=True)
# lattice.generate_hamiltonian(1, 2, 3)
