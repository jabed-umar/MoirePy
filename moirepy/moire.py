import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
import matplotlib.pyplot as plt

from .layers import Layer
from .utils import get_rotation_matrix, are_coeffs_integers
from . import moirepy_rust as rbck
from .utils import LatticeAlreadyFinalisedError

class BilayerMoireLattice:
    """Bilayer moire lattice wrapper with Rust-backed geometry and Hamiltonian assembly.

    This class builds two rotated/translating copies of a monolayer lattice,
    generates points in a common moire supercell, and exposes fast methods for
    connection generation, phase factors, and sparse Hamiltonian construction.

    Parameters
    ----------
    latticetype : type[Layer]
        Layer class (not instance) that subclasses :class:`moirepy.layers.Layer`.
    ll1, ll2 : int
        Integer coefficients defining the lower-layer supercell direction.
    ul1, ul2 : int
        Integer coefficients defining the upper-layer direction used to infer
        twist angle relative to the lower layer.
    n1, n2 : int, optional
        Number of supercell tiles along moire vectors ``mlv1`` and ``mlv2``.
    translate_upper : tuple[float, float], optional
        Translation applied to the upper layer after rotation.
    pbc : bool, optional
        Whether periodic boundary conditions are enabled.
    study_proximity : int, optional
        Neighbor-search shell depth used by each layer.
    k : int, optional
        Number of orbitals per site.
    verbose : bool, optional
        If ``True``, prints twist angle and generated-cell counts.

    Raises
    ------
    ValueError
        If ``latticetype`` is not a ``Layer`` subclass, or if the computed
        moire vectors are incompatible with layer periodicity.
    AssertionError
        If generated upper and lower point counts are inconsistent.
    """

    def __init__(
        self,
        latticetype: Layer,
        ll1: int, ll2: int,
        ul1: int, ul2: int,
        n1: int = 1, n2: int = 1,
        translate_upper=(0, 0),
        pbc: bool = True,
        study_proximity: int=1,
        k: int = 1,
        verbose=True,
    ):
        # make sure latticetype is a class, not an instance and class which inherits from Layer
        if not isinstance(latticetype, type) or not issubclass(latticetype, Layer):
            raise ValueError("latticetype must be a class that inherits from Layer.")

        # 1. Instantiate the Python-level Layers
        self.lower_lattice = latticetype(pbc=pbc, study_proximity=study_proximity)
        self.upper_lattice = latticetype(pbc=pbc, study_proximity=study_proximity)

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
        """Plot both layer point clouds and moire/simulation boundaries.

        Notes
        -----
        Adds artists to the current Matplotlib axes and sets equal aspect
        ratio. Call ``plt.show()`` or save the figure after invoking this
        method.
        """
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

    def generate_connections(self, inter_layer_radius: float = 1.0):
        """Populate connection instruction buffers for all hopping channels.

        Parameters
        ----------
        inter_layer_radius : float, optional
            Maximum distance for inter-layer neighbor search.

        Raises
        ------
        RuntimeError
            If lattice points have not been generated on either layer.
        """
        # 1. Ensure points exist on both layers before proceeding
        if self.lower_lattice.points is None or self.upper_lattice.points is None:
            raise RuntimeError("Lattice points must be generated before finding connections.")

        # 2. Call the Rust backend to fill HoppingInstruction buffers
        # This fills hop_ll, hop_uu, hop_ul, etc. in one high-speed pass.
        self._rust_class.generate_connections(float(inter_layer_radius))

    def get_phase(self, k):
        """Compute Peierls phase matrix for one or multiple wavevectors.

        Parameters
        ----------
        k : array-like
            Either shape ``(2,)`` for one wavevector or ``(N, 2)`` for a batch.

        Returns
        -------
        scipy.sparse.coo_matrix or list[scipy.sparse.coo_matrix]
            Phase matrix/matrices with shape ``(n_tot, n_tot)`` where
            ``n_tot = (n_lower + n_upper) * orbitals``.

        Raises
        ------
        ValueError
            If ``k`` does not have shape ``(2,)`` or ``(N, 2)``.
        """
        k = np.asarray(k, dtype=float)
        n_tot = (len(self.lower_lattice.points) + len(self.upper_lattice.points)) * self.orbitals

        if k.ndim == 1:
            if k.shape[0] != 2:
                raise ValueError(f"k must have shape (2,), got {k.shape}")
            data, rows, cols = self._rust_class.get_phase(float(k[0]), float(k[1]))
            return coo_matrix((data, (rows, cols)), shape=(n_tot, n_tot), dtype=np.complex128)

        if k.ndim == 2 and k.shape[1] == 2:
            if k.shape[0] == 0:
                raise ValueError("k with shape (N, 2) must have N > 0")

            out = []
            for i in range(k.shape[0]):
                data, rows, cols = self._rust_class.get_phase(float(k[i, 0]), float(k[i, 1]))
                out.append(coo_matrix((data, (rows, cols)), shape=(n_tot, n_tot), dtype=np.complex128))

            return out

        raise ValueError(f"k must have shape (2,) or (N, 2), got {k.shape}")

    def get_hopping_instructions(self):
        """Return raw hopping instruction objects from the Rust backend.

        Returns
        -------
        dict[str, object]
            Mapping with keys ``tll``, ``tuu``, ``tul``, ``tlu``,
            ``tlself``, and ``tuself``.
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

    def _prepare_hopping_input(self, val, instruction, layer_i, layer_j, extra_inputs=None):
        """Normalize pair-hopping input for Rust Hamiltonian builders.

        Parameters
        ----------
        val : float or callable
            Scalar hopping amplitude or callable
            ``val(coo_i, coo_j, R, labels_i, labels_j, moire, **extra_inputs)``
            returning an array with shape ``(N, k, k)``.
        instruction : object
            Rust instruction buffer with site/type indices.
        layer_i, layer_j : Layer
            Source and destination layers.
        extra_inputs : dict, optional
            Additional keyword arguments forwarded to ``val`` when callable.

        Returns
        -------
        float or np.ndarray
            Scalar value or flattened contiguous array of length
            ``N * k * k``.
        """
        if not callable(val):
            if val is None:
                return 0.0

            n_pairs = len(instruction.site_i)
            k = self.orbitals
            arr = np.asarray(val)

            # Scalar path: keep scalar so Rust uses Input::Scalar.
            if arr.ndim == 0:
                scalar = arr.item()
                return complex(scalar) if np.iscomplexobj(arr) else float(scalar)

            # Constant block path: (k, k) gets repeated for every pair.
            if arr.shape == (k, k):
                tiled = np.broadcast_to(arr, (n_pairs, k, k))
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(tiled, dtype=out_dtype).ravel()

            # Direct per-pair input path: (N, k, k)
            if arr.shape == (n_pairs, k, k):
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(arr, dtype=out_dtype).ravel()

            # Special-case k=1 convenience: accept shape (N,) as per-pair scalar vector.
            if k == 1 and arr.shape == (n_pairs,):
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(arr, dtype=out_dtype).ravel()

            raise ValueError(
                f"Invalid hopping input shape {arr.shape}; expected scalar, ({k},{k}), ({n_pairs},{k},{k})"
                + (f", or ({n_pairs},) for k=1" if k == 1 else "")
            )

        k = self.orbitals
        # 1. Map IDs to string labels
        basis_i = np.array(layer_i.basis_types)
        basis_j = np.array(layer_j.basis_types)
        labels_i = basis_i[np.array(instruction.type1)]
        labels_j = basis_j[np.array(instruction.type2)]

        # 2. Get coordinates
        coo_i = layer_i.points[instruction.site_i]
        coo_j = layer_j.points[instruction.site_j]

        # 3. Build R (supercell shift vectors), shape (N, 2)
        R = np.column_stack([instruction.dx, instruction.dy])

        # 4. Call user function. Expected output shape: (N, k, k)
        result = val(coo_i, coo_j, R, labels_i, labels_j, self, **extra_inputs if extra_inputs is not None else {})
        n_pairs = len(instruction.site_i)
        result_arr = np.asarray(result)

        # k=1 convenience: callable can return a scalar constant.
        if self.orbitals == 1 and result_arr.ndim == 0:
            scalar = result_arr.item()
            return complex(scalar) if np.iscomplexobj(result_arr) else float(scalar)

        # 5. Flatten to 1D array for Rust
        # Shape change: (N, k, k) -> (N * k * k,)
        out_dtype = np.complex128 if np.iscomplexobj(result_arr) else np.float64

        if result_arr.shape == (n_pairs, self.orbitals, self.orbitals):
            return np.ascontiguousarray(result_arr, dtype=out_dtype).ravel()

        if result_arr.shape == (self.orbitals, self.orbitals):
            tiled = np.broadcast_to(result_arr, (n_pairs, self.orbitals, self.orbitals))
            return np.ascontiguousarray(tiled, dtype=out_dtype).ravel()

        if self.orbitals == 1 and result_arr.shape == (n_pairs,):
            return np.ascontiguousarray(result_arr, dtype=out_dtype).ravel()

        raise ValueError(
            f"Callable hopping returned shape {result_arr.shape}; expected ({n_pairs},{self.orbitals},{self.orbitals})"
            + (f", ({self.orbitals},{self.orbitals}), scalar, or ({n_pairs},) for k=1" if self.orbitals == 1 else f" or ({self.orbitals},{self.orbitals})")
        )

    def _prepare_self_input(self, val, instruction, layer, extra_inputs=None):
        """Normalize on-site/self-channel input for Rust Hamiltonian builders.

        Parameters
        ----------
        val : float or callable
            Scalar onsite term or callable returning an array with shape
            ``(N, k, k)``.
        instruction : object
            Rust instruction buffer with site/type indices.
        layer : Layer
            Layer from which onsite terms are sampled.
        extra_inputs : dict, optional
            Additional keyword arguments forwarded to ``val`` when callable.

        Returns
        -------
        float or np.ndarray
            Scalar value or flattened contiguous array of length
            ``N * k * k``.
        """
        if not callable(val):
            if val is None:
                return 0.0

            n_sites = len(instruction.site_i)
            k = self.orbitals
            arr = np.asarray(val)

            # Scalar path: keep scalar so Rust uses Input::Scalar.
            if arr.ndim == 0:
                scalar = arr.item()
                return complex(scalar) if np.iscomplexobj(arr) else float(scalar)

            # Constant onsite block path: (k, k) repeated for every site.
            if arr.shape == (k, k):
                tiled = np.broadcast_to(arr, (n_sites, k, k))
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(tiled, dtype=out_dtype).ravel()

            # Direct per-site input path: (N, k, k)
            if arr.shape == (n_sites, k, k):
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(arr, dtype=out_dtype).ravel()

            # Special-case k=1 convenience: accept shape (N,) as per-site scalar vector.
            if k == 1 and arr.shape == (n_sites,):
                out_dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
                return np.ascontiguousarray(arr, dtype=out_dtype).ravel()

            raise ValueError(
                f"Invalid onsite input shape {arr.shape}; expected scalar, ({k},{k}), ({n_sites},{k},{k})"
                + (f", or ({n_sites},) for k=1" if k == 1 else "")
            )

        k = self.orbitals
        basis = np.array(layer.basis_types)
        labels_i = basis[np.array(instruction.ptype)]
        coo_i = layer.points[instruction.site_i]
        
        # Call user function. Expected output shape: (N, k, k)
        result = val(coo_i, labels_i, self, **extra_inputs if extra_inputs is not None else {})
        # result = np.array([
        #     val(c, l) for c, l in zip(coo_i, labels_i)
        # ])

        n_sites = len(instruction.site_i)
        result_arr = np.asarray(result)

        # k=1 convenience: callable can return a scalar constant.
        if self.orbitals == 1 and result_arr.ndim == 0:
            scalar = result_arr.item()
            return complex(scalar) if np.iscomplexobj(result_arr) else float(scalar)

        out_dtype = np.complex128 if np.iscomplexobj(result_arr) else np.float64

        if result_arr.shape == (n_sites, self.orbitals, self.orbitals):
            return np.ascontiguousarray(result_arr, dtype=out_dtype).ravel()

        if result_arr.shape == (self.orbitals, self.orbitals):
            tiled = np.broadcast_to(result_arr, (n_sites, self.orbitals, self.orbitals))
            return np.ascontiguousarray(tiled, dtype=out_dtype).ravel()

        if self.orbitals == 1 and result_arr.shape == (n_sites,):
            return np.ascontiguousarray(result_arr, dtype=out_dtype).ravel()

        raise ValueError(
            f"Callable onsite returned shape {result_arr.shape}; expected ({n_sites},{self.orbitals},{self.orbitals})"
            + (f", ({self.orbitals},{self.orbitals}), scalar, or ({n_sites},) for k=1" if self.orbitals == 1 else f" or ({self.orbitals},{self.orbitals})")
        )

    def generate_hamiltonian(
        self,
        tll=0.0, tuu=0.0, tlu=0.0, tul=0.0,
        tuself=0.0, tlself=0.0,
        data_type=np.complex128,
        extra_inputs=None
    ):
        """Assemble the tight-binding Hamiltonian as a sparse COO matrix.

        Parameters
        ----------
        tll, tuu, tlu, tul : float or callable, optional
            Intra-/inter-layer hopping terms. Each can be:
            (1) scalar, or (2) callable returning ``(N, k, k)`` blocks.
        tuself, tlself : float or callable, optional
            Upper/lower onsite terms. Same scalar-or-callable contract.
        data_type : np.dtype, optional
            Output sparse matrix dtype, typically ``np.float64`` or
            ``np.complex128``.
        extra_inputs : dict, optional
            Extra keyword arguments passed into user callables.

        Returns
        -------
        scipy.sparse.coo_matrix
            Hamiltonian matrix with shape ``(n_tot, n_tot)``.

        Notes
        -----
        Complex assembly is automatically selected if any input resolves to a
        complex array or if ``data_type`` is ``np.complex128``.
        """
        # 1. Get instructions from Rust
        instr = self.get_hopping_instructions()
        l_lat, u_lat = self.lower_lattice, self.upper_lattice

        # 2. Process all inputs into either floats or 1D arrays
        v_tlself = self._prepare_self_input(tlself, instr["tlself"], l_lat, extra_inputs)
        v_tuself = self._prepare_self_input(tuself, instr["tuself"], u_lat, extra_inputs)

        v_tll = self._prepare_hopping_input(tll, instr["tll"], l_lat, l_lat, extra_inputs)
        v_tuu = self._prepare_hopping_input(tuu, instr["tuu"], u_lat, u_lat, extra_inputs)
        v_tlu = self._prepare_hopping_input(tlu, instr["tlu"], l_lat, u_lat, extra_inputs)
        v_tul = self._prepare_hopping_input(tul, instr["tul"], u_lat, l_lat, extra_inputs)

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


    # --- Geometric Coefficients ---
    @property
    def ll1(self):
        """int: First lower-layer integer coefficient."""
        return self._rust_class.ll1
    @ll1.setter
    def ll1(self, value):
        raise LatticeAlreadyFinalisedError("ll1", self.__class__.__name__)

    @property
    def ll2(self):
        """int: Second lower-layer integer coefficient."""
        return self._rust_class.ll2
    @ll2.setter
    def ll2(self, value):
        raise LatticeAlreadyFinalisedError("ll2", self.__class__.__name__)

    @property
    def ul1(self):
        """int: First upper-layer integer coefficient."""
        return self._rust_class.ul1
    @ul1.setter
    def ul1(self, value):
        raise LatticeAlreadyFinalisedError("ul1", self.__class__.__name__)

    @property
    def ul2(self):
        """int: Second upper-layer integer coefficient."""
        return self._rust_class.ul2
    @ul2.setter
    def ul2(self, value):
        raise LatticeAlreadyFinalisedError("ul2", self.__class__.__name__)

    # --- Supercell Parameters ---
    @property
    def n1(self):
        """int: Number of moire tiles along :attr:`mlv1`."""
        return self._rust_class.n1
    @n1.setter
    def n1(self, value):
        raise LatticeAlreadyFinalisedError("n1", self.__class__.__name__)

    @property
    def n2(self):
        """int: Number of moire tiles along :attr:`mlv2`."""
        return self._rust_class.n2
    @n2.setter
    def n2(self, value):
        raise LatticeAlreadyFinalisedError("n2", self.__class__.__name__)

    @property
    def theta(self):
        """float: Twist angle in radians."""
        return self._rust_class.theta
    @theta.setter
    def theta(self, value):
        raise LatticeAlreadyFinalisedError("theta", self.__class__.__name__)

    # --- Global Simulation State ---
    @property
    def pbc(self):
        """bool: Whether periodic boundary conditions are enabled."""
        return self._rust_class.pbc
    @pbc.setter
    def pbc(self, value):
        raise LatticeAlreadyFinalisedError("pbc", self.__class__.__name__)

    @property
    def orbitals(self):
        """int: Number of orbitals per lattice site."""
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
