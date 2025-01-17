from typing import Union, Callable
import numpy as np
from layers import Layer, HexagonalLayer, SquareLayer, RhombusLayer, TriangularLayer, KagomeLayer
import matplotlib.pyplot as plt
from utils import get_rotation_matrix
import time

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

        plt.plot(*zip(*self.lower_lattice.points), 'ro')
        plt.plot(*zip(*self.upper_lattice.points), 'bo')

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

    def validate_input(self, a):
        if callable(a): return a
        return lambda theta: a

    def generate_hamiltonian(
        self,
        tuu: Union[float, int, Callable[[float], float]] = 1,
        tdd: Union[float, int, Callable[[float], float]] = 1,
        tdu: Union[float, int, Callable[[float], float]] = 1,
        k_upper_limit: int = 6,
    ):
        if isinstance(tuu, int) or isinstance(tuu, float): tuu = self.validate_input(tuu)
        if isinstance(tdd, int) or isinstance(tdd, float): tdd = self.validate_input(tdd)
        if isinstance(tdu, int) or isinstance(tdu, float): tdu = self.validate_input(tdu)
        assert callable(tuu) and callable(tdd) and callable(tdu), "tuu, tdd and tdu must be floats, ints or callable objects like functions"


        # self.plot_lattice()

        # 1. interaction inside the lower lattice
        ham_ll = np.zeros((len(self.lower_lattice.points), len(self.lower_lattice.points)))
        _, indices = self.lower_lattice.query_non_self(self.lower_lattice.points, k=6)
        for this_i in range(len(self.lower_lattice.points)):
            this_coo = self.lower_lattice.points[this_i]
            for neigh_i in indices[this_i]:
                neigh_coo = self.lower_lattice.points[neigh_i]
                # angle of this_coo to neigh_coo vector with respect to the x-axis
                angle = np.arctan2(neigh_coo[1] - this_coo[1], neigh_coo[0] - this_coo[0])
                ham_ll[this_i, neigh_i] = tuu(angle)
                # print(angle)
                # plot a line segment from this_coo to neigh_coo
        #         plt.plot([this_coo[0], neigh_coo[0]], [this_coo[1], neigh_coo[1]], 'k', linewidth=0.5)
        # plt.show()

        # 2. interaction inside the upper lattice
        ham_uu = np.zeros((len(self.upper_lattice.points), len(self.upper_lattice.points)))
        _, indices = self.upper_lattice.query_non_self(self.upper_lattice.points, k=6)
        for this_i in range(len(self.upper_lattice.points)):
            this_coo = self.upper_lattice.points[this_i]
            for neigh_i in indices[this_i]:
                neigh_coo = self.upper_lattice.points[neigh_i]
                # angle of this_coo to neigh_coo vector with respect to the x-axis
                angle = np.arctan2(neigh_coo[1] - this_coo[1], neigh_coo[0] - this_coo[0])
                ham_uu[this_i, neigh_i] = tdd(angle)
        #         plt.plot([this_coo[0], neigh_coo[0]], [this_coo[1], neigh_coo[1]], 'k', linewidth=0.5)
        # plt.show()

        # 3. interaction from the lower to the upper lattice
        ham_lu = np.zeros((len(self.lower_lattice.points), len(self.upper_lattice.points)))
        _, indices = self.lower_lattice.query(self.upper_lattice.points, k=1)
        for this_i in range(len(self.lower_lattice.points)):
            neigh_i = indices[this_i]
            ham_lu[this_i, neigh_i] = tdu(+1)

        # 4. interaction from the upper to the lower lattice
        ham_ul = np.zeros((len(self.upper_lattice.points), len(self.lower_lattice.points)))
        _, indices = self.upper_lattice.query(self.lower_lattice.points, k=1)
        for this_i in range(len(self.upper_lattice.points)):
            neigh_i = indices[this_i]
            ham_ul[this_i, neigh_i] = tdu(-1)

        # combine the hamiltonians
        self.ham = np.block([
            [ham_ll, ham_lu],
            [ham_ul, ham_uu]
        ])

        return self.ham


t = time.time()

# lattice = MoireLattice(TriangularLayer, 9, 10, 3+0, 2+0, pbc=False)
lattice = MoireLattice(TriangularLayer, 11, 12, 1, 1, pbc=True)
# lattice = MoireLattice(TriangularLayer, 12, 13, 1, 1, pbc=True)
# lattice = MoireLattice(TriangularLayer, 9, 10, 2, 4, pbc=False)

print(f"initialization took: {time.time() - t:.2f} seconds")
t = time.time()

ham = lattice.generate_hamiltonian(1, 1, 1)

print(f"hamiltonian generation took: {time.time() - t:.2f} seconds")

# check if ham is hermitian
if np.allclose(ham, ham.T.conj()): print("Hamiltonian is hermitian")
else: print("Hamiltonian is not hermitian")

t = time.time()

evals, evecs = np.linalg.eigh(ham)

print(f"diagonalization took: {time.time() - t:.2f} seconds")




# plt.imshow(ham)
# plt.colorbar()
# plt.show()

# lattice.plot_lattice()