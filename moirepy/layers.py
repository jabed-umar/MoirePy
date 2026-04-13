import numpy as np
import matplotlib.pyplot as plt
import importlib

from .utils import LatticeAlreadyFinalisedError


def _get_rust_backend():
    """Import the compiled Rust extension lazily to avoid package init cycles."""
    return importlib.import_module("moirepy.moirepy_rust")


class Layer:
    """Base class representing a 2D crystal lattice layer.

    Wraps a Rust-backed lattice object that handles the heavy numerical work.

    Parameters
    ----------
    lv1 : np.ndarray
        First primitive lattice vector. **Must** lie along the x-axis,
        i.e. ``lv1[1] == 0``.
    lv2 : np.ndarray
        Second primitive lattice vector. Must have a positive y-component
        (``lv2[1] > 0``).
    basis_points : list of tuple
        Atom positions within the unit cell. Each element is a 3-tuple
        ``(x, y, type_label)`` where ``x`` and ``y`` are fractional or
        Cartesian coordinates and ``type_label`` is a unique string
        identifier (e.g. ``"A"``, ``"B"``).
    neighbours : dict
        Nearest-neighbour displacement vectors for each atom type.
        Keys are the type labels from ``basis_points``; values are lists
        of ``[dx, dy]`` displacement vectors that reach the neighbours of
        that atom type.
    pbc : bool, optional
        Whether to apply periodic boundary conditions when computing
        neighbours. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells to include in the neighbour search.
        Larger values find neighbours across more unit-cell images.
        Default is ``1``.

    Raises
    ------
    ValueError
        If ``lv1`` is not aligned with the x-axis or ``lv2`` does not
        have a positive y-component.
    LatticeAlreadyFinalisedError
        If any of the lattice parameters are attempted to be reassigned after construction.

    See Also
    --------
    SquareLayer, TriangularLayer, HexagonalLayer, KagomeLayer, Rhombus60Layer :
        Pre-built layer types for common 2D lattices.

    Notes
    -----
    The orientation convention (``lv1`` along x, ``lv2`` with positive y)
    is required by the rotation machinery. See
    https://jabed-umar.github.io/MoirePy/find_theta/ for details.

    Examples
    --------
    Build a custom honeycomb-like layer from scratch:

    >>> import numpy as np
    >>> from moirepy import Layer
    >>> lv1 = np.array([1.0, 0.0])
    >>> lv2 = np.array([0.5, np.sqrt(3) / 2])
    >>> basis = [(0, 0, "A"), (1, 1/np.sqrt(3), "B")]
    >>> nb = {
    ...     "A": [[0, 1/np.sqrt(3)], [-0.5, -1/(2*np.sqrt(3))], [0.5, -1/(2*np.sqrt(3))]],
    ...     "B": [[0.5, 1/(2*np.sqrt(3))], [-0.5, 1/(2*np.sqrt(3))], [0, -1/np.sqrt(3)]],
    ... }
    >>> layer = Layer(lv1, lv2, basis, nb)
    """

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

        rbck = _get_rust_backend()
        self._rust_lattice = rbck.Layer(
            lv1=lv1.tolist(),
            lv2=lv2.tolist(),
            basis_data=[tuple(p) for p in basis_points],
            neighbours=neighbours,
            pbc=pbc,
            study_proximity=study_proximity
        )

    def perform_rotation_translation(self, rot: float, translation=None) -> None:
        """Apply a rotation and optional translation to the layer.

        Must be called **before** :meth:`generate_points`. The rotation is
        stored internally and applied when lattice points are generated.

        Parameters
        ----------
        rot : float
            Rotation angle in **degrees**.
        translation : array-like of length 2, optional
            Translation vector ``(tx, ty)`` applied after rotation.
            Defaults to ``(0.0, 0.0)`` if not provided.

        Raises
        ------
        RuntimeError
            If points have already been generated via :meth:`generate_points`.
        ValueError
            If ``translation`` cannot be interpreted as a pair of numbers.
        """
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
        """Populate the layer with lattice points inside a moiré supercell.

        The supercell is tiled ``mln1 x mln2`` times along ``mlv1`` and
        ``mlv2`` respectively. All basis atoms in each unit cell within that
        region are added to :attr:`points`.

        Parameters
        ----------
        mlv1 : np.ndarray
            First moiré supercell lattice vector, shape ``(2,)``.
        mlv2 : np.ndarray
            Second moiré supercell lattice vector, shape ``(2,)``.
        mln1 : int, optional
            Number of supercell repetitions along ``mlv1``. Must be a
            positive integer. Default is ``1``.
        mln2 : int, optional
            Number of supercell repetitions along ``mlv2``. Must be a
            positive integer. Default is ``1``.

        Raises
        ------
        AssertionError
            If ``mln1`` or ``mln2`` are not positive integers.
        """
        assert isinstance(mln1, int) and mln1 > 0, "mln1 must be a positive integer."
        assert isinstance(mln2, int) and mln2 > 0, "mln2 must be a positive integer."

        # Ensure inputs are numpy arrays with float type
        mlv1 = np.asarray(mlv1, dtype=float)
        mlv2 = np.asarray(mlv2, dtype=float)

        self._rust_lattice.generate_points(mlv1, mlv2, mln1, mln2)

    def first_nearest_neighbours(self, points=None, types=None):
        """Find the first nearest neighbours for each lattice point.

        For every point, returns the indices of its first-shell neighbours
        as defined by the ``neighbours`` displacement dictionary supplied
        at construction. The Rust backend handles type-to-ID mapping.

        Parameters
        ----------
        points : array-like of shape (N, 2), optional
            Query point coordinates. Defaults to :attr:`points` (the full
            generated lattice) if not provided.
        types : array-like of length N, optional
            Atom-type label for each query point. Defaults to
            :attr:`point_types` if not provided.

        Returns
        -------
        list of list of int
            A jagged list of length N. ``result[i]`` is the list of
            (global) indices of the first nearest neighbours of point ``i``.

        Raises
        ------
        ValueError
            If ``points`` is not a 2-D array with shape ``(N, 2)``.
        """
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
        """Find all lattice points within a given radius of each query point.

        Searches the full generated lattice for any atom within ``radius``
        of each entry in ``query_points``.

        Parameters
        ----------
        query_points : array-like of shape (N, 2)
            Coordinates of the N probe positions.
        radius : float
            Search radius (in the same units as the lattice vectors).

        Returns
        -------
        bigger_list : list of list of int
            ``bigger_list[i]`` contains the global indices of neighbours of
            query point ``i`` from the *larger* index perspective (used for
            deduplication in Hamiltonian assembly).
        smaller_list : list of list of int
            ``smaller_list[i]`` contains the corresponding partner indices
            (smaller-index end of each bond).
        dist_list : list of list of float
            ``dist_list[i]`` contains the Euclidean distances to each
            neighbour found for query point ``i``.

        Raises
        ------
        ValueError
            If ``query_points`` is not a 2-D array with shape ``(N, 2)``.

        Notes
        -----
        The ``bigger``/``smaller`` split mirrors the convention used in the
        Rust backend to avoid double-counting hopping terms when building
        tight-binding Hamiltonians.
        """
        query_points = np.asarray(query_points, dtype=float)
        if query_points.ndim != 2 or query_points.shape[1] != 2:
            raise ValueError(f"query_points must be (N, 2), got {query_points.shape}")

        n_query = len(query_points)
        q_indices, bigger_indices, smaller_indices, distances = \
            self._rust_lattice.get_neighbors_within_radius(query_points, float(radius))

        # Group flat arrays into per-query-point lists of lists.
        # q_indices is like [0, 0, 0, 1, 1, 2, 2, 2, ...] — one entry per match.
        bigger_list  = [[] for _ in range(n_query)]
        smaller_list = [[] for _ in range(n_query)]
        dist_list    = [[] for _ in range(n_query)]
        for qi, bi, si, d in zip(q_indices, bigger_indices, smaller_indices, distances):
            bigger_list[qi].append(int(bi))
            smaller_list[qi].append(int(si))
            dist_list[qi].append(float(d))

        return bigger_list, smaller_list, dist_list

    def plot_lattice(self, plot_connections: bool = True, colours: list = ["r", "g", "b", "c", "m", "y", "k"]) -> None:
        """Plot the generated lattice points, optionally with neighbour bonds.

        Atom types are distinguished by colour, cycling through the
        ``colours`` list in the same order as :attr:`basis_points`. Bond
        lines are drawn with 10 % opacity so they do not obscure the sites.

        Parameters
        ----------
        plot_connections : bool, optional
            Whether to draw lines between each atom and its nearest
            neighbours. Default is ``True``.
        colours : list of str, optional
            Matplotlib colour strings assigned to each atom type in order.
            If the list has exactly one entry, all types share that colour.
            Default is ``["r", "g", "b", "c", "m", "y", "k"]``.

        Notes
        -----
        This method adds artists to the *current* active Matplotlib axes.
        Call ``plt.show()`` (or save the figure) after this method.
        The axes aspect ratio is set to ``"equal"`` automatically.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> layer = HexagonalLayer()
        >>> layer.generate_points(mlv1, mlv2)
        >>> layer.plot_lattice()
        >>> plt.show()
        """
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
        """np.ndarray: First primitive lattice vector ``[x, 0]``."""
        return self._rust_lattice.lv1

    @lv1.setter
    def lv1(self, value):
        raise LatticeAlreadyFinalisedError("lv1", self.__class__.__name__)

    @property
    def lv2(self):
        """np.ndarray: Second primitive lattice vector with positive y-component."""
        return self._rust_lattice.lv2

    @lv2.setter
    def lv2(self, value):
        raise LatticeAlreadyFinalisedError("lv2", self.__class__.__name__)

    @property
    def basis_points(self):
        """list of tuple: Atom positions ``(x, y, type_label)`` within the unit cell."""
        return self._rust_lattice.basis_points

    @basis_points.setter
    def basis_points(self, value):
        raise LatticeAlreadyFinalisedError("basis_points", self.__class__.__name__)

    @property
    def basis_types(self):
        """list of str: Ordered atom-type labels extracted from ``basis_points``."""
        return self._rust_lattice.basis_types

    @basis_types.setter
    def basis_types(self, value):
        raise LatticeAlreadyFinalisedError("basis_types", self.__class__.__name__)

    @property
    def neighbours(self):
        """dict: Nearest-neighbour displacements keyed by atom-type label.

        Returns a plain Python dict mapping each type label (str) to its
        list of ``[dx, dy]`` displacement vectors. Reconstructed from the
        Rust backend on every access.
        """
        return dict(zip(self.basis_types, self._rust_lattice.neighbours))

    @neighbours.setter
    def neighbours(self, value):
        raise LatticeAlreadyFinalisedError("neighbours", self.__class__.__name__)

    @property
    def mlv1(self):
        """np.ndarray or None: First moiré supercell vector, set by :meth:`generate_points`."""
        return self._rust_lattice.mlv1

    @mlv1.setter
    def mlv1(self, value):
        raise LatticeAlreadyFinalisedError("mlv1", self.__class__.__name__)

    @property
    def mlv2(self):
        """np.ndarray or None: Second moiré supercell vector, set by :meth:`generate_points`."""
        return self._rust_lattice.mlv2

    @mlv2.setter
    def mlv2(self, value):
        raise LatticeAlreadyFinalisedError("mlv2", self.__class__.__name__)

    @property
    def mln1(self):
        """int: Supercell repetition count along ``mlv1``."""
        return self._rust_lattice.mln1

    @mln1.setter
    def mln1(self, value):
        raise LatticeAlreadyFinalisedError("mln1", self.__class__.__name__)

    @property
    def mln2(self):
        """int: Supercell repetition count along ``mlv2``."""
        return self._rust_lattice.mln2

    @mln2.setter
    def mln2(self, value):
        raise LatticeAlreadyFinalisedError("mln2", self.__class__.__name__)

    @property
    def pbc(self):
        """bool: Whether periodic boundary conditions are active."""
        return self._rust_lattice.pbc

    @pbc.setter
    def pbc(self, value):
        raise LatticeAlreadyFinalisedError("pbc", self.__class__.__name__)

    @property
    def study_proximity(self):
        """int: Number of supercell shells included in the neighbour search."""
        return self._rust_lattice.study_proximity

    @study_proximity.setter
    def study_proximity(self, value):
        raise LatticeAlreadyFinalisedError("study_proximity", self.__class__.__name__)

    @property
    def points(self):
        """np.ndarray or None: Generated lattice-point coordinates, shape ``(N, 2)``.

        ``None`` until :meth:`generate_points` has been called.
        """
        return self._rust_lattice.points

    @points.setter
    def points(self, value):
        raise LatticeAlreadyFinalisedError("points", self.__class__.__name__)

    @property
    def bigger_points(self):
        """np.ndarray or None: Lattice points that extend beyond the primary supercell.

        Used internally to resolve periodic images during neighbour search.
        ``None`` until :meth:`generate_points` has been called.
        """
        return self._rust_lattice.bigger_points

    @bigger_points.setter
    def bigger_points(self, value):
        raise LatticeAlreadyFinalisedError("bigger_points", self.__class__.__name__)

    @property
    def point_types(self):
        """list of str: Atom-type label for each point in :attr:`points`.

        Has the same length as :attr:`points`. ``None`` until
        :meth:`generate_points` has been called.
        """
        return self._rust_lattice.point_types

    @point_types.setter
    def point_types(self, value):
        raise LatticeAlreadyFinalisedError("point_types", self.__class__.__name__)

    @property
    def rot_m(self):
        """np.ndarray: 2x2 rotation matrix applied to this layer."""
        return self._rust_lattice.rot_m

    @rot_m.setter
    def rot_m(self, value):
        raise LatticeAlreadyFinalisedError("rot_m", self.__class__.__name__)

    @property
    def translation(self):
        """tuple of float: Translation ``(tx, ty)`` applied after rotation."""
        return self._rust_lattice.translation

    @translation.setter
    def translation(self, value):
        raise LatticeAlreadyFinalisedError("translation", self.__class__.__name__)


# ===============================================
# ============= some example layers =============
# ===============================================

class SquareLayer(Layer):
    """A simple square lattice with one atom per unit cell.

    Lattice vectors are ``[1, 0]`` and ``[0, 1]``. The single atom type
    ``"A"`` has four nearest neighbours at the cardinal directions.

    Parameters
    ----------
    pbc : bool, optional
        Periodic boundary conditions. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells for the neighbour search. Default is ``1``.

    Examples
    --------
    >>> from moirepy import SquareLayer
    >>> layer = SquareLayer()
    """

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
    """A triangular (hexagonal close-packed) lattice with one atom per unit cell.

    Lattice vectors are ``[1, 0]`` and ``[0.5, sqrt(3)/2]``. The single atom
    type ``"A"`` has six nearest neighbours at the vertices of a regular
    hexagon.

    Parameters
    ----------
    pbc : bool, optional
        Periodic boundary conditions. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells for the neighbour search. Default is ``1``.

    Examples
    --------
    >>> from moirepy import TriangularLayer
    >>> layer = TriangularLayer()
    """

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
    """A rhombic lattice with a 60-degree opening angle and one atom per cell.

    Both lattice vectors have unit length; ``lv1`` is along x and ``lv2``
    is at 60 degrees. The atom type ``"A"`` has four nearest neighbours
    (two along each lattice direction).

    Parameters
    ----------
    pbc : bool, optional
        Periodic boundary conditions. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells for the neighbour search. Default is ``1``.

    Notes
    -----
    The opening angle is hard-coded to 60 degrees. For other angles, derive
    a subclass or use :class:`Layer` directly.

    Examples
    --------
    >>> from moirepy import Rhombus60Layer
    >>> layer = Rhombus60Layer()
    """

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
    """A Kagome lattice with three atoms per unit cell (types ``"A"``, ``"B"``, ``"C"``).

    The Kagome lattice is a corner-sharing triangular tiling. Each atom type
    has four nearest neighbours. Lattice vectors are ``[1, 0]`` and
    ``[0.5, sqrt(3)/2]``.

    Parameters
    ----------
    pbc : bool, optional
        Periodic boundary conditions. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells for the neighbour search. Default is ``1``.

    Examples
    --------
    >>> from moirepy import KagomeLayer
    >>> layer = KagomeLayer()
    """

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
    """A honeycomb (hexagonal) lattice with two atoms per unit cell (types ``"A"`` and ``"B"``).

    This is the lattice structure of graphene and related 2D materials.
    Each atom has three nearest neighbours of the opposite type. Lattice
    vectors are ``[1, 0]`` and ``[0.5, sqrt(3)/2]``.

    Parameters
    ----------
    pbc : bool, optional
        Periodic boundary conditions. Default is ``False``.
    study_proximity : int, optional
        Number of supercell shells for the neighbour search. Default is ``1``.

    Examples
    --------
    >>> from moirepy import HexagonalLayer
    >>> layer = HexagonalLayer()
    """

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
