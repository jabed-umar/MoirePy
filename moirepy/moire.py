from typing import Union, Callable
import numpy as np
from layers import Layer, SquareLayer, Rhombus60Layer, TriangularLayer
import matplotlib.pyplot as plt
from utils import get_rotation_matrix
import time

class MoireLattice:
    def __init__(self, latticetype: Layer, a:int, b:int, nx:int=1, ny:int=1, interlayer_distance=1, pbc = False):
        # study_proximity = 1 means only studying nearest neighbours will be eabled, 2 means study of next nearest neighbours will be enabled too and so on
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
        self.interlayer_distance = interlayer_distance

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

    def _validate_input(self, a, name):
        if a is None:
            a = 1
            print(f"WARNING: {name} is not set, setting it to 1")
        if callable(a): return a
        return lambda this_coo, neigh_coo: a
    

    def generate_hamiltonian(
        self,
        tuu: Union[float, int, Callable[[float], float]] = None,
        tll: Union[float, int, Callable[[float], float]] = None,
        tlu: Union[float, int, Callable[[float], float]] = None,
        tul: Union[float, int, Callable[[float], float]] = None
    ):
        if tuu is None or isinstance(tuu, int) or isinstance(tuu, float): tuu = self._validate_input(tuu, "tuu")
        if tll is None or isinstance(tll, int) or isinstance(tll, float): tll = self._validate_input(tll, "tll")
        if tlu is None or isinstance(tlu, int) or isinstance(tlu, float): tlu = self._validate_input(tlu, "tlu")
        if tul is None or isinstance(tul, int) or isinstance(tul, float): tul = self._validate_input(tul, "tul")
        assert callable(tuu) and callable(tll) and callable(tlu) and callable(tul), "tuu, tll, tlu and tul must be floats, ints or callable objects like functions"
        # self.plot_lattice()

        # 1. interaction inside the lower lattice
        ham_ll = np.zeros((len(self.lower_lattice.points), len(self.lower_lattice.points)))
        _, indices = self.lower_lattice.first_nearest_neighbours(self.lower_lattice.points, self.lower_lattice.point_types)
        for this_i in range(len(self.lower_lattice.points)):
            this_coo = self.lower_lattice.points[this_i]
            for neigh_i in indices[this_i]:
                neigh_coo = self.lower_lattice.points[neigh_i]
                ham_ll[this_i, neigh_i] = tuu(this_coo, neigh_coo)

        # 2. interaction inside the upper lattice
        ham_uu = np.zeros((len(self.upper_lattice.points), len(self.upper_lattice.points)))
        _, indices = self.upper_lattice.first_nearest_neighbours(self.upper_lattice.points, self.upper_lattice.point_types)
        for this_i in range(len(self.upper_lattice.points)):
            this_coo = self.upper_lattice.points[this_i]
            for neigh_i in indices[this_i]:
                neigh_coo = self.upper_lattice.points[neigh_i]
                ham_uu[this_i, neigh_i] = tll(this_coo, neigh_coo)

        # 3. interaction from the lower to the upper lattice
        ham_lu = np.zeros((len(self.lower_lattice.points), len(self.upper_lattice.points)))
        _, indices = self.upper_lattice.query(self.lower_lattice.points, k=1)
        for this_i in range(len(self.lower_lattice.points)):
            neigh_i = indices[this_i]
            this_coo = self.lower_lattice.points[this_i]
            neigh_coo = self.upper_lattice.points[neigh_i]
            ham_lu[this_i, neigh_i] = tlu(this_coo, neigh_coo)

        # 4. interaction from the upper to the lower lattice
        ham_ul = np.zeros((len(self.upper_lattice.points), len(self.lower_lattice.points)))
        _, indices = self.lower_lattice.query(self.upper_lattice.points, k=1)
        for this_i in range(len(self.upper_lattice.points)):
            neigh_i = indices[this_i]
            this_coo = self.upper_lattice.points[this_i]
            neigh_coo = self.lower_lattice.points[neigh_i]
            ham_ul[this_i, neigh_i] = tul(this_coo, neigh_coo)

        # in ham_ll and ham_uu, check if sum of all the rows is same
        print(f"unique sums in ham_ll: {np.unique(np.sum(ham_ll, axis=1))}")
        print(f"unique sums in ham_uu: {np.unique(np.sum(ham_uu, axis=1))}")
        print(f"unique sums in ham_lu: {np.unique(np.sum(ham_lu, axis=1))}")
        print(f"unique sums in ham_ul: {np.unique(np.sum(ham_ul, axis=1))}")

        # combine the hamiltonians
        self.ham = np.block([
            [ham_ll, ham_lu],
            [ham_ul, ham_uu]
        ])


        return self.ham


t = time.time()

# lattice = MoireLattice(TriangularLayer, 9, 10, 3+0, 2+0, pbc=False)
# lattice = MoireLattice(TriangularLayer, 5, 6, 2, 2, pbc=True)
lattice = MoireLattice(SquareLayer, 6, 7, 2, 2, pbc=True)
# lattice = MoireLattice(TriangularLayer, 5, 6, 2, 2, pbc=True)
# lattice = MoireLattice(TriangularLayer, 12, 13, 1, 1, pbc=True)
# lattice = MoireLattice(TriangularLayer, 9, 10, 2, 4, pbc=False)

print(f"initialization took: {time.time() - t:.2f} seconds")
t = time.time()

ham = lattice.generate_hamiltonian(1, 1, 1)

print(f"hamiltonian generation took: {time.time() - t:.2f} seconds")

# check if ham is hermitian
if np.allclose(ham, ham.T.conj()): print("Hamiltonian is hermitian.")
else: print("Hamiltonian is not hermitian.")

t = time.time()

evals, evecs = np.linalg.eigh(ham)

print(f"diagonalization took: {time.time() - t:.2f} seconds")




plt.imshow(ham, cmap="gray")
plt.colorbar()
plt.show()

# lattice.plot_lattice()