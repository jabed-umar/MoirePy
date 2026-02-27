from typing import Union, Callable
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
import matplotlib.pyplot as plt

from .layers import Layer
from .utils import get_rotation_matrix, are_coeffs_integers
from . import moirepy_rust as rbck
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
        # self.lower_lattice.plot_lattice(colours=["b"], plot_connections=True)
        # self.upper_lattice.plot_lattice(colours=["r"], plot_connections=True)

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

    def _prepare_hopping_input(self, val, instruction, layer_i, layer_j):
        """Prepares tll, tuu, tlu, tul. Handles scalars or (N, k, k) blocks."""
        if not callable(val):
            # If scalar, we still send a float. Rust will handle the 'Scalar' variant.
            return float(val) if val is not None else 0.0

        k = self.orbitals
        # 1. Map IDs to string labels
        basis_i = np.array(layer_i.basis_types)
        basis_j = np.array(layer_j.basis_types)
        labels_i = basis_i[np.array(instruction.type1)]
        labels_j = basis_j[np.array(instruction.type2)]

        # 2. Get coordinates
        coo_i = layer_i.points[instruction.site_i]
        coo_j = layer_j.points[instruction.site_j]

        # 3. Call user function. Expected output shape: (N, k, k)
        # result = val(coo_i, coo_j, labels_i, labels_j)
        result = np.array([
            val(ci, cj, li, lj) for ci, cj, li, lj in zip(coo_i, coo_j, labels_i, labels_j)
        ])

        # 4. Flatten to 1D array for Rust
        # Shape change: (N, k, k) -> (N * k * k,)
        return np.ascontiguousarray(result, dtype=np.complex128 if np.iscomplexobj(result) else np.float64).flatten()

    def _prepare_self_input(self, val, instruction, layer):
        """Prepares tlself, tuself. Handles scalars or (N, k, k) blocks."""
        if not callable(val):
            return float(val) if val is not None else 0.0

        k = self.orbitals
        basis = np.array(layer.basis_types)
        labels_i = basis[np.array(instruction.ptype)]
        coo_i = layer.points[instruction.site_i]
        
        # Call user function. Expected output shape: (N, k, k)
        # result = val(coo_i, labels_i)
        result = np.array([
            val(c, l) for c, l in zip(coo_i, labels_i)
        ])

        return np.ascontiguousarray(result, dtype=np.complex128 if np.iscomplexobj(result) else np.float64).flatten()

    def generate_hamiltonian(
        self,
        tll=0.0, tuu=0.0, tlu=0.0, tul=0.0,
        tuself=0.0, tlself=0.0,
        data_type=np.float64
    ):
        # 1. Get instructions from Rust
        instr = self.get_hopping_instructions()
        l_lat, u_lat = self.lower_lattice, self.upper_lattice

        # 2. Process all inputs into either floats or 1D arrays
        v_tlself = self._prepare_self_input(tlself, instr["tlself"], l_lat)
        v_tuself = self._prepare_self_input(tuself, instr["tuself"], u_lat)

        v_tll = self._prepare_hopping_input(tll, instr["tll"], l_lat, l_lat)
        v_tuu = self._prepare_hopping_input(tuu, instr["tuu"], u_lat, u_lat)
        v_tlu = self._prepare_hopping_input(tlu, instr["tlu"], l_lat, u_lat)
        v_tul = self._prepare_hopping_input(tul, instr["tul"], u_lat, l_lat)

        # 3. Determine if we need the complex or real Rust gatekeeper
        # If any prepared input is an array of complex numbers, we must use complex
        is_complex = any(np.iscomplexobj(v) for v in (v_tll, v_tuu, v_tlu, v_tul, v_tuself, v_tlself))
        is_complex = is_complex or (data_type == np.complex128)

        if is_complex:
            data, rows, cols = self._rust_class.build_ham_complex(
                v_tll, v_tuu, v_tlu, v_tul, v_tuself, v_tlself
            )
        else:
            data, rows, cols = self._rust_class.build_ham(
                v_tll, v_tuu, v_tlu, v_tul, v_tuself, v_tlself
            )

        # 4. Construct the matrix
        n_tot = (len(l_lat.points) + len(u_lat.points)) * self.orbitals
        return coo_matrix((data, (rows, cols)), shape=(n_tot, n_tot), dtype=data_type)




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
        suppress_obc_warning: bool = False,
    ):
        if not suppress_nxny_warning and (self.n1 != 1 or self.n2 != 1):
            raise ValueError("WARNING: n1 or n2 != 1. Momentum space is usually for n1=n2=1."
                "Aborting mission. set suppress_nxny_warning=True to override this check.")

        if not suppress_obc_warning and not self.pbc:
            raise ValueError(
                "k-space Hamiltonian requested for an OBC system. "
                "Since all lattice shifts R=0, the result is identical to the real-space matrix and k-independent. "
                "Set suppress_obc_warning=True to proceed anyway."
            )

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

