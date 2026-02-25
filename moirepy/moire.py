from typing import Union, Callable, Sequence
import numpy as np
from .layers import Layer
import matplotlib.pyplot as plt
from .utils import get_rotation_matrix, are_coeffs_integers
from . import moirepy_rust as rbck  # import rust backend
from scipy.sparse import coo_matrix
from .utils import get_rotation_matrix, LatticeAlreadyFinalisedError

class COOBuilder:
    def __init__(self, rows=None, cols=None, data=None):
        self.rows = rows if rows is not None else []
        self.cols = cols if cols is not None else []
        self.data = data if data is not None else []
        assert len(self.rows) == len(self.cols) == len(self.data), "Initial rows, cols, and data must be of the same length."

    def add(self, r, c, val):
        self.rows.append(r)
        self.cols.append(c)
        self.data.append(val)
        



from time import perf_counter_ns

class Timer:
    def __init__(self):
        self.last = perf_counter_ns()
        self.times = []
    def lap(self):
        now = perf_counter_ns()
        d = now - self.last
        self.times.append(d/1e6)
        self.last = now



class BilayerMoireLattice:
    def __init__(
        self,
        latticetype: Layer,
        ll1: int, ll2: int,
        ul1: int, ul2: int,
        n1: int = 1, n2: int = 1,
        translate_upper=(0, 0),
        pbc: bool = True,
        k: int = 1,
        verbose=True,
    ):
        # make sure latticetype is a class, not an instance and class which inherits from Layer
        if not isinstance(latticetype, type) or not issubclass(latticetype, Layer):
            raise ValueError("latticetype must be a class that inherits from Layer.")
        
        # 1. Instantiate the Python-level Layers
        self.lower_lattice = latticetype(pbc=pbc)
        self.upper_lattice = latticetype(pbc=pbc)

        # 2. Local variables for geometric setup
        lv1, lv2 = self.lower_lattice.lv1, self.lower_lattice.lv2
        c = np.dot(lv1, lv2) / (np.linalg.norm(lv1) * np.linalg.norm(lv2))
        beta = np.arccos(c)
        
        mlv1 = ll1 * lv1 + ll2 * lv2
        mlv2 = get_rotation_matrix(beta).dot(mlv1)

        one = ll1 * lv1 + ll2 * lv2
        two = ul1 * lv1 + ul2 * lv2
        
        # Calculate twist angle
        c_theta = np.dot(one, two) / (np.linalg.norm(one) * np.linalg.norm(two))
        theta = np.arccos(np.clip(c_theta, -1.0, 1.0))
        
        if verbose: 
            print(f"twist angle = {theta:.4f} rad ({np.rad2deg(theta):.4f} deg)")

        # 3. Synchronize Layer state
        self.upper_lattice.perform_rotation_translation(theta, translate_upper)
        
        # Validation
        if not (are_coeffs_integers(self.lower_lattice.lv1, self.lower_lattice.lv2, mlv1) and
                are_coeffs_integers(self.upper_lattice.lv1, self.upper_lattice.lv2, mlv1)):
            raise ValueError("Moiré lattice vectors are inconsistent with layer periodicity.")

        # Trigger point generation (The 22ms native loops)
        self.lower_lattice.generate_points(mlv1, mlv2, n1, n2)
        self.upper_lattice.generate_points(mlv1, mlv2, n1, n2)

        # 4. Instantiate and store the Rust backend

        self._rust_class = rbck.BilayerMoire(
            self.lower_lattice._rust_lattice, # Passing Rust handles directly
            self.upper_lattice._rust_lattice,
            ll1, ll2, ul1, ul2,
            n1, n2,
            theta,
            mlv1, mlv2,
            np.array(translate_upper, dtype=np.float64),
            pbc,
            k
        )

        if verbose:
            print(f"{len(self.upper_lattice.points)} cells in upper lattice")
            print(f"{len(self.lower_lattice.points)} cells in lower lattice")
        assert len(self.lower_lattice.points) == len(self.upper_lattice.points), "FATAL ERROR: number of cells in lower and upper lattice are not equal, report and try different ll1, ll2, ul1, ul2 values."

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

    def generate_connections(self, inter_layer_radius: float = 3.0):
        """
        Populates the internal Rust hopping buffers for intra- and inter-layer connections.
        
        Args:
            inter_layer_radius (float): The cutoff distance for hopping between layers.
        """
        # 1. Ensure points exist on both layers before proceeding
        if self.lower_lattice.points is None or self.upper_lattice.points is None:
            raise RuntimeError("Lattice points must be generated before finding connections.")

        # 2. Call the Rust backend to fill HoppingInstruction buffers
        # This fills hop_ll, hop_uu, hop_ul, etc. in one high-speed pass.
        self._rust_class.generate_connections(float(inter_layer_radius))

    def get_hopping_instructions(self):
        """
        Retrieves the raw connection data from the Rust backend.
        Returns a dictionary containing site indices and type IDs for all hopping terms.
        """
        # Assuming we expose these via getters in moire.rs
        return {
            "tll": self._rust_class.hop_ll,
            "tuu": self._rust_class.hop_uu,
            "tul": self._rust_class.hop_ul,
            "tlu": self._rust_class.hop_lu,
            "tlself": self._rust_class.hop_l,
            "tuself": self._rust_class.hop_u,
        }



    def generate_hamiltonian(
        self,
        tll: float = 0.0,
        tuu: float = 0.0,
        tlu: float = 0.0,
        tul: float = 0.0,
        tuself: float = 0.0,
        tlself: float = 0.0,
        inter_layer_radius: float = 3.0,
        data_type=np.float64,
    ):
        t = Timer()
        
        # 1. Ensure all inputs are floats for the Rust backend
        tll, tuu, tlu, tul, tuself, tlself = [
            float(t) if t is not None else 0.0 
            for t in (tll, tuu, tlu, tul, tuself, tlself)
        ]
        
        data, rows, cols = self._rust_class.build_ham_from_scalars(
            tll, tuu, tlu, tul, tuself, tlself
        )
        
        t.lap()  # 1

        n_lower = len(self.lower_lattice.points)
        n_upper = len(self.upper_lattice.points)
        total_dim = (n_lower + n_upper) * self.orbitals
        
        ham = coo_matrix(
            (data, (rows, cols)),
            shape=(total_dim, total_dim),
            dtype=data_type
        )
        
        t.lap()  # 2
        
        return ham, t.times





    # def _as_callable(self, val, n_args=4):
    #     if callable(val):
    #         return val
    #     if n_args == 4:
    #         return lambda c1, c2, t1, t2: val if val is not None else 0
    #     return lambda c, t: val if val is not None else 0
    
    # def _validate_hamiltonian_inputs(self, tll, tuu, tlu, tul, tlself, tuself):
    #     """Helper to ensure all hopping terms are callable."""
    #     tll, tuu, tlu, tul = [self._as_callable(t, 4) for t in (tll, tuu, tlu, tul)]
    #     tlself, tuself = [self._as_callable(t, 2) for t in (tlself, tuself)]
    #     assert all(
    #         callable(fn)
    #         for fn in (tll, tuu, tlu, tul, tlself, tuself)
    #     ), "Hopping parameters must be floats, ints, or callable functions."
    #     return tll, tuu, tlu, tul, tlself, tuself

    # def generate_hamiltonian(
    #     self,
    #     tll: Union[float, int, Callable] = None,
    #     tuu: Union[float, int, Callable] = None,
    #     tlu: Union[float, int, Callable] = None,
    #     tul: Union[float, int, Callable] = None,
    #     tlself: Union[float, int, Callable] = None,
    #     tuself: Union[float, int, Callable] = None,
    #     inter_layer_radius: float = 3.0,
    #     data_type: np.dtype = np.float64,
    # ):
    #     k = self.orbitals
    #     n_lower = len(self.lower_lattice.points)
    #     n_upper = len(self.upper_lattice.points)
    #     total_dim = (n_lower + n_upper) * k
        
    #     tll, tuu, tlu, tul, tlself, tuself = self._validate_hamiltonian_inputs(tll, tuu, tlu, tul, tlself, tuself)
    #     builder = COOBuilder()

    #     # 1. Lower Lattice Intra-layer (Indices: 0 to n_lower*k - 1)
    #     for i in range(n_lower):
    #         val = tlself(self.lower_lattice.points[i], self.lower_lattice.point_types[i])
    #         for o in range(k):
    #             builder.add(i*k + o, i*k + o, val)

    #     bigger_indices, indices, _ = self.lower_lattice.first_nearest_neighbours(self.lower_lattice.points, self.lower_lattice.point_types)
    #     for this_i in range(n_lower):
    #         this_coo, this_type = self.lower_lattice.points[this_i], self.lower_lattice.point_types[this_i]
    #         for phantom_neigh_i, neigh_i in zip(bigger_indices[this_i], indices[this_i]):
    #             neigh_coo = self.lower_lattice.bigger_points[phantom_neigh_i] if self.pbc else self.lower_lattice.points[neigh_i]
    #             neigh_type = self.lower_lattice.point_types[neigh_i]
    #             val = tll(this_coo, neigh_coo, this_type, neigh_type)
    #             for o1 in range(k):
    #                 for o2 in range(k):
    #                     builder.add(this_i*k + o1, neigh_i*k + o2, val)

    #     # 2. Upper Lattice Intra-layer (Indices: n_lower*k to total_dim - 1)
    #     offset = n_lower * k
    #     for i in range(n_upper):
    #         val = tuself(self.upper_lattice.points[i], self.upper_lattice.point_types[i])
    #         for o in range(k):
    #             builder.add(offset + i*k + o, offset + i*k + o, val)
            
    #     bigger_indices, indices, _ = self.upper_lattice.first_nearest_neighbours(self.upper_lattice.points, self.upper_lattice.point_types)
    #     for this_i in range(n_upper):
    #         this_coo, this_type = self.upper_lattice.points[this_i], self.upper_lattice.point_types[this_i]
    #         for phantom_neigh_i, neigh_i in zip(bigger_indices[this_i], indices[this_i]):
    #             neigh_coo = self.upper_lattice.bigger_points[phantom_neigh_i] if self.pbc else self.upper_lattice.points[neigh_i]
    #             neigh_type = self.upper_lattice.point_types[neigh_i]
    #             val = tuu(this_coo, neigh_coo, this_type, neigh_type)
    #             for o1 in range(k):
    #                 for o2 in range(k):
    #                     builder.add(offset + this_i*k + o1, offset + neigh_i*k + o2, val)

    #     # 3. Inter-layer (Radius Search)
    #     u_indices, l_indices, l_coords = self.lower_lattice.get_neighbors_within_radius(self.upper_lattice.points, inter_layer_radius)
        
    #     for i in range(len(u_indices)):
    #         u_idx, l_idx = u_indices[i], l_indices[i]
    #         u_pos, l_pos = self.upper_lattice.points[u_idx], l_coords[i]
    #         u_type, l_type = self.upper_lattice.point_types[u_idx], self.lower_lattice.point_types[l_idx]
            
    #         # tul: Upper -> Lower hopping
    #         val_ul = tul(u_pos, l_pos, u_type, l_type)
    #         # tlu: Lower -> Upper hopping (assuming Hermiticity if not provided)
    #         val_lu = tlu(l_pos, u_pos, l_type, u_type)

    #         for o1 in range(k):
    #             for o2 in range(k):
    #                 # Upper row, Lower col (H_ul block)
    #                 builder.add(offset + u_idx*k + o1, l_idx*k + o2, val_ul)
    #                 # Lower row, Upper col (H_lu block)
    #                 builder.add(l_idx*k + o1, offset + u_idx*k + o2, val_lu)

        # from scipy.sparse import coo_matrix
        # self.ham = coo_matrix((builder.data, (builder.rows, builder.cols)), shape=(total_dim, total_dim), dtype=data_type)
        # return self.ham.tocsc()

    def generate_k_space_hamiltonian(
        self,
        k: np.ndarray,
        tll: Union[float, int, Callable] = None,
        tuu: Union[float, int, Callable] = None,
        tlu: Union[float, int, Callable] = None,
        tul: Union[float, int, Callable] = None,
        tlself: Union[float, int, Callable] = None,
        tuself: Union[float, int, Callable] = None,
        inter_layer_radius: float = 3.0,
        suppress_nxny_warning: bool = False,
        suppress_pbc_warning: bool = False,
    ):
        if not suppress_nxny_warning and (self.n1 != 1 or self.n2 != 1):
            print("WARNING: n1 or n2 != 1. Momentum space is usually for n1=n2=1. Aborting mission. set suppress_nxny_warning=True to override this check.")
            return

        if not suppress_pbc_warning and not self.pbc:
            print("WARNING: k-space generation is only physically meaningful with pbc=True. Aborting mission. set suppress_pbc_warning=True to override this check.")
            return

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


    # --- Geometric Coefficients ---
    @property
    def ll1(self):
        return self._rust_class.ll1
    @ll1.setter
    def ll1(self, value):
        raise LatticeAlreadyFinalisedError("ll1", self.__class__.__name__)

    @property
    def ll2(self):
        return self._rust_class.ll2
    @ll2.setter
    def ll2(self, value):
        raise LatticeAlreadyFinalisedError("ll2", self.__class__.__name__)

    @property
    def ul1(self):
        return self._rust_class.ul1
    @ul1.setter
    def ul1(self, value):
        raise LatticeAlreadyFinalisedError("ul1", self.__class__.__name__)

    @property
    def ul2(self):
        return self._rust_class.ul2
    @ul2.setter
    def ul2(self, value):
        raise LatticeAlreadyFinalisedError("ul2", self.__class__.__name__)

    # --- Supercell Parameters ---
    @property
    def n1(self):
        return self._rust_class.n1
    @n1.setter
    def n1(self, value):
        raise LatticeAlreadyFinalisedError("n1", self.__class__.__name__)

    @property
    def n2(self):
        return self._rust_class.n2
    @n2.setter
    def n2(self, value):
        raise LatticeAlreadyFinalisedError("n2", self.__class__.__name__)

    @property
    def theta(self):
        return self._rust_class.theta
    @theta.setter
    def theta(self, value):
        raise LatticeAlreadyFinalisedError("theta", self.__class__.__name__)

    # --- Global Simulation State ---
    @property
    def pbc(self):
        return self._rust_class.pbc
    @pbc.setter
    def pbc(self, value):
        raise LatticeAlreadyFinalisedError("pbc", self.__class__.__name__)

    @property
    def orbitals(self):
        return self._rust_class.orbitals
    @orbitals.setter
    def orbitals(self, value):
        raise LatticeAlreadyFinalisedError("orbitals", self.__class__.__name__)

    @property
    def mlv1(self):
        """First Moiré lattice vector."""
        return self._rust_class.mlv1
    @mlv1.setter
    def mlv1(self, value):
        raise LatticeAlreadyFinalisedError("mlv1", self.__class__.__name__)

    @property
    def mlv2(self):
        """Second Moiré lattice vector."""
        return self._rust_class.mlv2
    @mlv2.setter
    def mlv2(self, value):
        raise LatticeAlreadyFinalisedError("mlv2", self.__class__.__name__)

    @property
    def translate_upper(self):
        """Translation vector applied to the upper layer."""
        return self._rust_class.translate_upper
    @translate_upper.setter
    def translate_upper(self, value):
        raise LatticeAlreadyFinalisedError("translate_upper", self.__class__.__name__)


# class BilayerMPMoireLattice(BilayerMoireLattice):
#     def __init__(
#         self,
#         latticetype: Layer,
#         ll1:int, ll2:int,  # lower lattice
#         ul1:int, ul2:int,  # upper lattice
#         n1:int=1, n2:int=1,
#         translate_upper=(0, 0),
#         pbc:bool=True,
#         k:int=1,  # number of orbitals
#     ):
#         pass

