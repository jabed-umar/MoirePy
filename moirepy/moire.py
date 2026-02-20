from typing import Union, Callable, Sequence
import numpy as np
from .layers import Layer
import matplotlib.pyplot as plt
from .utils import get_rotation_matrix, are_coeffs_integers

class COOBuilder:
    def __init__(self, rows=[], cols=[], data=[]):
        self.rows, self.cols, self.data = rows, cols, data
        assert len(self.rows) == len(self.cols) == len(self.data), "Initial rows, cols, and data must be of the same length."

    def add(self, r, c, val):
        self.rows.append(r)
        self.cols.append(c)
        self.data.append(val)

class BilayerMoireLattice:  # both layers same, only one point in one unit cell
    def __init__(
        self,
        latticetype: Layer,
        ll1:int, ll2:int,  # lower lattice
        ul1:int, ul2:int,  # upper lattice
        n1:int=1, n2:int=1,
        translate_upper=(0, 0),
        pbc:bool=True,
        k:int=1,  # number of orbitals
        verbose = True,
    ):
        """
        Initializes a Moiré lattice composed of two twisted layers of the same type.

        Args:
            latticetype (Layer): A subclass of the `Layer` class representing the lattice type used for both layers.
            ll1, ll2, ul1, ul2 (int): Values select from the [AVC tool](https://jabed-umar.github.io/MoirePy/theory/avc/).
            n1 (int, optional): Number of moiré cells along the first lattice vector.
            n2 (int, optional): Number of moiré cells along the second lattice vector.
            translate_upper (tuple, optional): Translation vector (dx, dy) applied to the upper layer before rotation.
            pbc (bool, optional): Whether to apply periodic boundary conditions. If False, open boundary conditions are used.
            k (int, optional): Number of orbitals on each lattice point.
        """

        # study_proximity = 1 means only studying nearest neighbours will be enabled,
        # 2 means study of next nearest neighbours will be enabled too and so on,
        # always better to keep this value 1 or two more than what you will actually need.
        lower_lattice = latticetype(pbc=pbc)
        upper_lattice = latticetype(pbc=pbc)

        lv1, lv2 = lower_lattice.lv1, lower_lattice.lv2

        # c = cos(theta) between lv1 and lv2 (60 degree for triangular, 90 for square and so on)
        c = np.dot(lv1, lv2) / (np.linalg.norm(lv1) * np.linalg.norm(lv2))
        beta = np.arccos(c)
        mlv1 = ll1*lv1 + ll2*lv2  # because lower latice is fixed
        mlv2 = get_rotation_matrix(beta).dot(mlv1)

        # calculating the moire twist angle
        one = ll1*lv1 + ll2*lv2  # the coords of overlapping point in the lower lattice
        two = ul1*lv1 + ul2*lv2  # the coords of overlapping point in the upper lattice
        assert np.isclose(np.linalg.norm(one), np.linalg.norm(two)), "INPUT ERROR: the two points are not overlapping, check ll1, ll2, ul1, ul2 values"
        c = np.dot(one, two) / (np.linalg.norm(one) * np.linalg.norm(two))
        theta = np.arccos(c)  # in radians
        if verbose: print(f"twist angle = {theta:.4f} rad ({np.rad2deg(theta):.4f} deg)")

        upper_lattice.perform_rotation_translation(theta, translate_upper)
        assert (
            are_coeffs_integers(lower_lattice.lv1, lower_lattice.lv2, mlv1) and
            are_coeffs_integers(upper_lattice.lv1, upper_lattice.lv2, mlv1)
        ), "FATAL ERROR: calculated mlv2 is incorrect"
        lower_lattice.generate_points(mlv1, mlv2, n1, n2)
        upper_lattice.generate_points(mlv1, mlv2, n1, n2)
        # print(f"{mlv1 = }")
        # print(f"{mlv2 = }")

        self.ll1 = ll1
        self.ll2 = ll2
        self.ul1 = ul1
        self.ul2 = ul2
        self.n1 = n1
        self.n2 = n2
        self.translate_upper = translate_upper
        self.lower_lattice = lower_lattice
        self.upper_lattice = upper_lattice
        self.theta = theta
        self.mlv1 = mlv1
        self.mlv2 = mlv2
        self.pbc = pbc
        self.orbitals = k
        self.ham = None

        if verbose:
            print(f"{len(self.upper_lattice.points)} cells in upper lattice")
            print(f"{len(self.lower_lattice.points)} cells in lower lattice")
        assert len(self.lower_lattice.points) == len(self.upper_lattice.points), "FATAL ERROR: number of cells in lower and upper lattice are not equal, report and take different ll1, ll2, ul1, ul2 values"

        # self.plot_lattice()

    def plot_lattice(self):
        mlv1 = self.mlv1
        mlv2 = self.mlv2
        n1 = self.n1
        n2 = self.n2

        plt.plot(*zip(*self.lower_lattice.points), 'r.', markersize=2)
        plt.plot(*zip(*self.upper_lattice.points), 'b.', markersize=2)
        self.lower_lattice.plot_lattice(colours=["b"], plot_connections=True)
        self.upper_lattice.plot_lattice(colours=["r"], plot_connections=True)

        # parallellogram around the whole lattice
        plt.plot([0, n1*mlv1[0]], [0, n1*mlv1[1]], 'k', linewidth=1)
        plt.plot([0, n2*mlv2[0]], [0, n2*mlv2[1]], 'k', linewidth=1)
        plt.plot([n1*mlv1[0], n1*mlv1[0] + n2*mlv2[0]], [n1*mlv1[1], n1*mlv1[1] + n2*mlv2[1]], 'k', linewidth=1)
        plt.plot([n2*mlv2[0], n1*mlv1[0] + n2*mlv2[0]], [n2*mlv2[1], n1*mlv1[1] + n2*mlv2[1]], 'k', linewidth=1)

        # just plot mlv1 and mlv2 parallellogram
        plt.plot([0, mlv1[0]], [0, mlv1[1]], 'k', linewidth=1)
        plt.plot([0, mlv2[0]], [0, mlv2[1]], 'k', linewidth=1)
        plt.plot([mlv1[0], mlv1[0] + mlv2[0]], [mlv1[1], mlv1[1] + mlv2[1]], 'k', linewidth=1)
        plt.plot([mlv2[0], mlv1[0] + mlv2[0]], [mlv2[1], mlv1[1] + mlv2[1]], 'k', linewidth=1)

        # set equal aspect ratio
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.grid()
        # plt.show()
        # plt.savefig("moire.pdf", bbox_inches='tight')

    def _as_callable(self, val, n_args=4):
        if callable(val):
            return val
        if n_args == 4:
            return lambda c1, c2, t1, t2: val if val is not None else 0
        return lambda c, t: val if val is not None else 0
    
    def _validate_hamiltonian_inputs(self, tll, tuu, tlu, tul, tuself, tlself):
        """Helper to ensure all hopping terms are callable."""
        tll, tuu, tlu, tul = [self._as_callable(t, 4) for t in (tll, tuu, tlu, tul)]
        tuself, tlself = [self._as_callable(t, 2) for t in (tuself, tlself)]
        assert all(
            callable(fn)
            for fn in (tll, tuu, tlu, tul, tuself, tlself)
        ), "Hopping parameters must be floats, ints, or callable functions."
        return tll, tuu, tlu, tul, tuself, tlself

    def generate_hamiltonian(
        self,
        tll: Union[float, int, Callable] = None,
        tuu: Union[float, int, Callable] = None,
        tlu: Union[float, int, Callable] = None,
        tul: Union[float, int, Callable] = None,
        tuself: Union[float, int, Callable] = None,
        tlself: Union[float, int, Callable] = None,
        inter_layer_radius: float = 3.0,
        data_type: np.dtype = np.float64,
    ):
        k = self.orbitals
        n_lower = len(self.lower_lattice.points)
        n_upper = len(self.upper_lattice.points)
        total_dim = (n_lower + n_upper) * k
        
        tll, tuu, tlu, tul, tuself, tlself = self._validate_hamiltonian_inputs(tll, tuu, tlu, tul, tuself, tlself)
        builder = COOBuilder()

        # 1. Lower Lattice Intra-layer (Indices: 0 to n_lower*k - 1)
        for i in range(n_lower):
            val = tlself(self.lower_lattice.points[i], self.lower_lattice.point_types[i])
            for o in range(k):
                builder.add(i*k + o, i*k + o, val)

        bigger_indices, indices, _ = self.lower_lattice.first_nearest_neighbours(self.lower_lattice.points, self.lower_lattice.point_types)
        for this_i in range(n_lower):
            this_coo, this_type = self.lower_lattice.points[this_i], self.lower_lattice.point_types[this_i]
            for phantom_neigh_i, neigh_i in zip(bigger_indices[this_i], indices[this_i]):
                neigh_coo = self.lower_lattice.bigger_points[phantom_neigh_i] if self.pbc else self.lower_lattice.points[neigh_i]
                neigh_type = self.lower_lattice.point_types[neigh_i]
                val = tll(this_coo, neigh_coo, this_type, neigh_type)
                for o1 in range(k):
                    for o2 in range(k):
                        builder.add(this_i*k + o1, neigh_i*k + o2, val)

        # 2. Upper Lattice Intra-layer (Indices: n_lower*k to total_dim - 1)
        offset = n_lower * k
        for i in range(n_upper):
            val = tuself(self.upper_lattice.points[i], self.upper_lattice.point_types[i])
            for o in range(k):
                builder.add(offset + i*k + o, offset + i*k + o, val)
            
        bigger_indices, indices, _ = self.upper_lattice.first_nearest_neighbours(self.upper_lattice.points, self.upper_lattice.point_types)
        for this_i in range(n_upper):
            this_coo, this_type = self.upper_lattice.points[this_i], self.upper_lattice.point_types[this_i]
            for phantom_neigh_i, neigh_i in zip(bigger_indices[this_i], indices[this_i]):
                neigh_coo = self.upper_lattice.bigger_points[phantom_neigh_i] if self.pbc else self.upper_lattice.points[neigh_i]
                neigh_type = self.upper_lattice.point_types[neigh_i]
                val = tuu(this_coo, neigh_coo, this_type, neigh_type)
                for o1 in range(k):
                    for o2 in range(k):
                        builder.add(offset + this_i*k + o1, offset + neigh_i*k + o2, val)

        # 3. Inter-layer (Radius Search)
        u_indices, l_indices, l_coords = self.lower_lattice.get_neighbors_within_radius(self.upper_lattice.points, inter_layer_radius)
        
        for i in range(len(u_indices)):
            u_idx, l_idx = u_indices[i], l_indices[i]
            u_pos, l_pos = self.upper_lattice.points[u_idx], l_coords[i]
            u_type, l_type = self.upper_lattice.point_types[u_idx], self.lower_lattice.point_types[l_idx]
            
            # tul: Upper -> Lower hopping
            val_ul = tul(u_pos, l_pos, u_type, l_type)
            # tlu: Lower -> Upper hopping (assuming Hermiticity if not provided)
            val_lu = tlu(l_pos, u_pos, l_type, u_type)

            for o1 in range(k):
                for o2 in range(k):
                    # Upper row, Lower col (H_ul block)
                    builder.add(offset + u_idx*k + o1, l_idx*k + o2, val_ul)
                    # Lower row, Upper col (H_lu block)
                    builder.add(l_idx*k + o1, offset + u_idx*k + o2, val_lu)

        from scipy.sparse import coo_matrix
        self.ham = coo_matrix((builder.data, (builder.rows, builder.cols)), shape=(total_dim, total_dim), dtype=data_type)
        return self.ham.tocsc()

    def generate_k_space_hamiltonian(
        self,
        k: np.ndarray,
        tll: Union[float, int, Callable] = None,
        tuu: Union[float, int, Callable] = None,
        tlu: Union[float, int, Callable] = None,
        tul: Union[float, int, Callable] = None,
        tuself: Union[float, int, Callable] = None,
        tlself: Union[float, int, Callable] = None,
        inter_layer_radius: float = 3.0,
        suppress_nxny_warning: bool = False,
    ):
        if not suppress_nxny_warning and (self.n1 != 1 or self.n2 != 1):
            print("WARNING: n1 or n2 != 1. Momentum space is usually for n1=n2=1.")

        # Validate and convert inputs to callables
        tll, tuu, tlu, tul, tuself, tlself = self._validate_hamiltonian_inputs(
            tll, tuu, tlu, tul, tuself, tlself
        )

        # Bloch phase factor helper
        def part(this_coo, neigh_coo):
            return np.exp(1j * (k @ (this_coo.squeeze() - neigh_coo.squeeze())))

        return self.generate_hamiltonian(
            tll=lambda c1, c2, t1, t2: tll(c1, c2, t1, t2) * part(c1, c2),
            tuu=lambda c1, c2, t1, t2: tuu(c1, c2, t1, t2) * part(c1, c2),
            tlu=lambda c1, c2, t1, t2: tlu(c1, c2, t1, t2) * part(c1, c2),
            tul=lambda c1, c2, t1, t2: tul(c1, c2, t1, t2) * part(c1, c2),
            tuself=tuself,
            tlself=tlself,
            inter_layer_radius=inter_layer_radius,
            data_type=np.complex128
        )



class BilayerMPMoireLattice(BilayerMoireLattice):
    def __init__(
        self,
        latticetype: Layer,
        ll1:int, ll2:int,  # lower lattice
        ul1:int, ul2:int,  # upper lattice
        n1:int=1, n2:int=1,
        translate_upper=(0, 0),
        pbc:bool=True,
        k:int=1,  # number of orbitals
    ):
        pass

