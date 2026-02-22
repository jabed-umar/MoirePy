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
        lattice_points: list,
        neighbours: dict,
        pbc: bool = False,
        study_proximity: int = 1
    ) -> None:
        self._rust_lattice = rbck.Layer(
            lv1=lv1.tolist(),
            lv2=lv2.tolist(),
            basis_data=[tuple(p) for p in lattice_points],
            neighbours=neighbours,
            pbc=pbc,
            study_proximity=study_proximity
        )

    def perform_rotation_translation(self, rot: float, translation: np.ndarray = None) -> None:
        # Check against the Rust backend directly
        if self._rust_lattice.points.shape[0] > 0:
            raise RuntimeError("Cannot perform rotation and translation after points have been generated.")
        
        # Handle default value and force type conversion
        if translation is None:
            translation = np.array([0.0, 0.0])
        translation = np.asarray(translation, dtype=float)

        # Call the Rust backend
        self._rust_lattice.perform_rotation_translation(rot, translation)


    # def generate_points(
    #         self,
    #         mlv1: np.ndarray,
    #         mlv2: np.ndarray,
    #         mln1: int = 1,
    #         mln2: int = 1,
    # ) -> None:
    #     """
    #     Generates points for a Moiré lattice based on the given lattice
    #     vectors and the number of unit cells along each direction.

    #     Args:
    #         mlv1 (np.ndarray): The first Moiré lattice vector.
    #         mlv2 (np.ndarray): The second Moiré lattice vector.
    #         mln1 (int, optional): The number of Moiré unit cells
    #             along the first lattice vector. Defaults to 1.
    #         mln2 (int, optional): The number of Moiré unit
    #             cells along the second lattice vector. Defaults to 1.

    #     Returns:
    #         None: The function modifies the object state and
    #         stores the generated points and their types.

    #     Raises:
    #         AssersionError: If mln1 and mln2 are not positive integers.
    #         AssersionError: if test is True and mln1 or mln2 is not 1

    #     Example:
    #     ```python
    #     >>> lattice = MyLattice()  # a class inheriting the Layer class
    #     >>> lattice.generate_points(np.array([1, 0]), np.array([0.5, np.sqrt(3)/2]), mln1=1, mln2=1)
    #     >>> print(lattice.points)
    #     ```
    #     """
        
    #     # raise the promised errors:
    #     assert isinstance(mln1, int) and mln1 > 0, "mln1 must be a positive integer."
    #     assert isinstance(mln2, int) and mln2 > 0, "mln2 must be a positive integer."
        
    #     self.mlv1 = mlv1  # Moire lattice vector 1
    #     self.mlv2 = mlv2  # Moire lattice vector 2
    #     self.mln1 = mln1  # Number of moire unit cells along mlv1
    #     self.mln2 = mln2  # Number of moire unit cells along mlv2

    #     # Step 1: Find the maximum distance to determine the number of points along each direction
    #     points = [np.array([0, 0]), mlv1, mlv2, mlv1 + mlv2]
    #     max_distance = max(
    #         np.linalg.norm(points[0] - points[1]),
    #         np.linalg.norm(points[0] - points[2]),
    #         np.linalg.norm(points[0] - points[3]),
    #     )

    #     # Calculate number of grid points based on maximum distance and lattice vectors
    #     n = math.ceil(max_distance / min(np.linalg.norm(self.lv1), np.linalg.norm(self.lv2))) * 2

    #     # Step 2: Generate points inside one moire unit cell (based on `lv1` and `lv2`)
    #     step1_points = []  # List to hold points inside the unit cell
    #     step1_names = []  # List to hold the names of the points
    #     for i in range(-n, n + 1):  # Iterate along mlv1
    #         for j in range(-n, n + 1):  # Iterate along mlv2
    #             # Calculate the lattice point inside the unit cell
    #             point_o = i * self.lv1 + j * self.lv2
    #             for xpos, ypos, name in self.lattice_points:
    #                 point = point_o + np.array([xpos, ypos])
    #                 step1_points.append(point)
    #                 step1_names.append(name)

    #     step1_points = np.array(step1_points)
    #     step1_names = np.array(step1_names)

    #     # Apply the boundary check method (inside_boundaries) to filter the points
    #     mask = self._inside_boundaries(step1_points, 1, 1)  # generating points only inside the first cell now
    #     step1_points = step1_points[mask]
    #     step1_names = step1_names[mask]

    #     # Step 3: Copy and translate the unit cell to create the full lattice
    #     points = []  # List to hold all the moire points
    #     names = []
    #     for i in range(self.mln1):  # Translate along mlv1 direction
    #         for j in range(self.mln2):  # Translate along mlv2 direction
    #             translation_vector = i * mlv1 + j * mlv2
    #             translated_points = step1_points + translation_vector  # Translate points
    #             points.append(translated_points)
    #             names.append(step1_names)

    #     self.points = np.vstack(points)
    #     self.point_types = np.hstack(names)
    #     self.generate_kdtree()

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

    # def _point_positions(self, points: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    #     """
    #     Determines the position of each point relative to a parallelogram defined by vectors A and B.

    #     Args:
    #         points (np.ndarray): Array of points to be analyzed.
    #         A (np.ndarray): The first vector of the parallelogram.
    #         B (np.ndarray): The second vector of the parallelogram.

    #     Returns:
    #         np.ndarray: An array indicating the position of each point:
    #             - (0, 0) for points inside the parallelogram.
    #             - (-1, 1) or (1, -1) for points outside on specific sides.
    #             left side and right side will give -1 and 1 respectively
    #             top side and bottom side will give -1 and 1 respectively
    #     """
    #     # Compute determinants for positions relative to OA and BC
    #     det_OA = (points[:, 0] * A[1] - points[:, 1] * A[0]) <= self.toll_scale * 1e-2
    #     det_BC = ((points[:, 0] - B[0]) * A[1] - (points[:, 1] - B[1]) * A[0]) <= self.toll_scale * 1e-2
    #     position_y = det_OA.astype(float) + det_BC.astype(float)

    #     # Compute determinants for positions relative to OB and AC
    #     det_OB = (points[:, 0] * B[1] - points[:, 1] * B[0]) > -self.toll_scale * 1e-2
    #     det_AC = ((points[:, 0] - A[0]) * B[1] - (points[:, 1] - A[1]) * B[0]) > -self.toll_scale * 1e-2
    #     position_x = det_OB.astype(float) + det_AC.astype(float)

    #     return np.column_stack((position_x, position_y)) - 1

    # def _inside_polygon(self, points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    #     """
    #     Determines if each point is inside a polygon using the ray-casting method.


    #     Args:
    #         points (np.ndarray): Array of points to check.
    #         polygon (np.ndarray): Vertices of the polygon to
    #             test against, in counterclockwise order.

    #     Returns:
    #         np.ndarray: A boolean array where True indicates
    #             that the point is inside the polygon.

    #     Example:
    #     ```python
    #     points = np.array([[0.5, 0.5], [1, 1], [-1, -1]])
    #     polygon = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    #     inside = _inside_polygon(points, polygon)
    #     print(inside)
    #     ```
    #     """
    #     x, y = points[:, 0], points[:, 1]
    #     px, py = polygon[:, 0], polygon[:, 1]
    #     px_next, py_next = np.roll(px, -1), np.roll(py, -1)
    #     edge_cond = (y[:, None] > np.minimum(py, py_next)) & (y[:, None] <= np.maximum(py, py_next))
    #     with np.errstate(divide='ignore', invalid='ignore'):
    #         xinters = np.where(py != py_next, (y[:, None] - py) * (px_next - px) / (py_next - py) + px, np.inf)
    #     ray_crosses = edge_cond & (x[:, None] <= xinters)
    #     inside = np.sum(ray_crosses, axis=1) % 2 == 1
    #     return inside  # mask

    # def _inside_boundaries(self, points: np.ndarray, mln1=None, mln2=None) -> np.ndarray:
    #     """
    #     Determines if the given points lie within the boundaries of the Moiré lattice pattern.

    #     Args:
    #         points (np.ndarray): Array of points to check.
    #         mln1 (int, optional): The number of unit cells along the first direction.
    #             Defaults to the object's current value.
    #         mln2 (int, optional): The number of unit cells along the second direction.
    #             Defaults to the object's current value.

    #     Returns:
    #         np.ndarray: A boolean array where True indicates that
    #             the point is within the boundaries of the lattice.

    #     Raises:
    #         ValueError: If the points array has an invalid shape.

    #     Example:
    #     ```python
    #     points = np.array([[0.5, 0.5], [2, 2], [-1, -1]])
    #     lattice_boundaries = _inside_boundaries(points, mln1=3, mln2=3)
    #     print(lattice_boundaries)
    #     ```
    #     """
    #     v1 = (mln1 if mln1 else self.mln1) * self.mlv1
    #     v2 = (mln2 if mln2 else self.mln2) * self.mlv2

    #     p1 = np.array([0, 0])
    #     p2 = np.array([v1[0], v1[1]])
    #     p3 = np.array([v2[0], v2[1]])
    #     p4 = np.array([v1[0] + v2[0], v1[1] + v2[1]])

    #     shift_dir = -(v1 + v2)
    #     shift_dir = shift_dir / np.linalg.norm(shift_dir)  # normalize
    #     shift = shift_dir * self.toll_scale * 1e-4

    #     return self._inside_polygon(
    #         points,
    #         np.array([p1, p2, p4, p3]) + shift
    #     )

    # def generate_kdtree(self) -> None:
    #     """
    #     Generates a KDTree for spatial queries of points in the Moiré lattice.
    #     If PBC is enabled, additional points outside the primary unit cell are
    #     considered for accurate queries (same numbers of neigbours for all points).

    #     Returns:
    #         None: The function modifies the object state by generating
    #             a KDTree for spatial queries.

    #     Raises:
    #         ValueError: If the points in the lattice are not defined.
    #     """
    #     if not self.pbc:  # OBC is easy
    #         self.kdtree = KDTree(self.points)
    #         return

    #     # in case of periodic boundary conditions, we need to generate a bigger set of points
    #     all_points = []
    #     all_point_names = []
    #     for i in range(-1, 2):
    #         for j in range(-1, 2):
    #             all_points.append(self.points + i * self.mln1 * self.mlv1 + j * self.mln2 * self.mlv2)
    #             all_point_names.append(self.point_types)

    #     all_points = np.vstack(all_points)
    #     all_point_names = np.hstack(all_point_names)

    #     v1 = self.mln1 * self.mlv1
    #     v2 = self.mln2 * self.mlv2

    #     neigh_pad_1 = (1 + self.study_proximity) * np.linalg.norm(self.lv1) / np.linalg.norm(v1)
    #     neigh_pad_2 = (1 + self.study_proximity) * np.linalg.norm(self.lv2) / np.linalg.norm(v2)

    #     mask = self._inside_polygon(all_points, np.array([
    #         (-neigh_pad_1) * v1 + (-neigh_pad_2) * v2,
    #         (1 + neigh_pad_1) * v1 + (-neigh_pad_2) * v2,
    #         (1 + neigh_pad_1) * v1 + (1 + neigh_pad_2) * v2,
    #         (-neigh_pad_1) * v1 + (1 + neigh_pad_2) * v2,
    #     ]))
    #     # print(mask.shape, mask.dtype)
    #     points = all_points[mask]
    #     point_names = all_point_names[mask]
    #     self.bigger_points = points
    #     self.bigger_point_types = point_names
    #     self.kdtree = KDTree(points)
    #     self._generate_mapping()

    # def _generate_mapping(self) -> None:
    #     """
    #     Generates a mapping between two sets of points (larger and smaller lattices) based on their positions
    #     and computes the distances between corresponding points. If the distance between a point in the larger lattice
    #     and its nearest neighbor in the smaller lattice exceeds a specified tolerance, it raises a `ValueError` and plots
    #     the lattice points for visualization.

    #     This function uses a KDTree to find the nearest neighbor in the smaller lattice for each point in the larger lattice.
    #     It stores the resulting mappings in the `self.mappings` dictionary, where keys are indices in `self.bigger_points`
    #     and values are the corresponding indices in `self.points`.

    #     Raises:
    #         ValueError: If the distance between a point and its nearest neighbor exceeds the tolerance defined by `self.toll_scale`.

    #     Example:
    #     ```python
    #     my_lattice._generate_mapping()
    #     ```
    #     The function performs the following steps:
    #     1. Initializes an empty dictionary `self.mappings`.
    #     2. Uses a KDTree to query the neighbors for each point in the larger lattice (`self.bigger_points`).
    #     3. Computes the translation needed for each point based on a lattice scaling factor.
    #     4. If the distance between the corresponding points exceeds the tolerance, it raises a `ValueError` and plots the points.
    #     5. Stores the index mappings of the larger lattice points to smaller lattice points in `self.mappings`.
    #     6. The lattice plots show the parallelograms formed by `mlv1` and `mlv2` vectors for visualization.

    #     Note:
    #         - The function assumes `self.points` and `self.bigger_points`
    #         are defined as numpy arrays with the coordinates of the points in the lattices.
    #         - The translations used in the function are calculated based on `mln1`, `mln2`,
    #         `mlv1`, and `mlv2`, which define the lattice scaling and vectors.

    #     """
    #     self.mappings = {}
    #     tree = KDTree(self.points)  # smaller set
    #     translations = self._point_positions(
    #         self.bigger_points,
    #         self.mln1 * self.mlv1,
    #         self.mln2 * self.mlv2
    #     )
    #     for i, (dx, dy) in enumerate(translations):
    #         mapped_point = self.bigger_points[i] - (dx * self.mlv1 * self.mln1 + dy * self.mlv2 * self.mln2)
    #         distance, index = tree.query(mapped_point)
    #         if distance >= self.toll_scale * 1e-3:

    #             raise ValueError(f"FATAL ERROR: Distance {distance} exceeds tolerance for {i}th point {self.bigger_points[i]} mapped at location {mapped_point} with translation ({dx}, {dy}).")
    #         self.mappings[i] = index

    #     # point positions... for each point in self.point, point position is a array of length 2 (x, y)
    #     # where the elemnts are -1, 0 and 1... this is what their value mean about their position
    #     # (-1, 1) | (0, 1) | (1, 1)
    #     # -----------------------------
    #     # (-1, 0) | (0, 0) | (1, 0)
    #     # -----------------------------
    #     # (-1,-1) | (0,-1) | (1,-1)
    #     # -----------------------------
    #     # (0, 0) is our actual lattice part...
    #     # do this for all points in self.bigger_points:
    #     # all point with point_positions = (x, y) need to be translated by
    #     # (-x*self.mlv1*self.mln1 - y*self.mlv2*self.mln2) to get the corresponding point inside the lattice
    #     # then you would need to run a query on a newly kdtree of the smaller points...
    #     # to the get the index of the corresponding point inside the lattice (distance should be zero, just saying)
    #     # now we already know the index of the point in the self.bigger_points... so we can map that to the index of the point in the self.points
    #     # then we will store that in `self.mappings``
    #     # self.mapppings will be a dictionary with keys as the indices in the
    #     # self.bigger_points (unique) and values as the indices in the self.points (not unique)


    def first_nearest_neighbours(self, points: np.ndarray, types: np.ndarray):
        assert self.kdtree is not None
        
        unique_types = np.unique(types)
        max_neighs = max(len(self.neighbours[t]) for t in unique_types)
        n_points = points.shape[0]

        bigger_indices_arr = np.full((n_points, max_neighs), -1, dtype=int)
        smaller_indices_arr = np.full((n_points, max_neighs), -1, dtype=int)
        smaller_distances_arr = np.full((n_points, max_neighs), np.nan, dtype=float)

        # Pre-compute mapping array for vectorized indexing
        if self.pbc:
            map_arr = np.array([self.mappings[i] for i in range(len(self.bigger_points))])

        for t in unique_types:
            # 1. Mask to find all points of this specific type
            mask = (types == t)
            type_points = points[mask]  # Shape (N_t, 2)
            rel_neighs = np.array(self.neighbours[t])  # Shape (M_t, 2)
            n_curr = len(rel_neighs)

            # 2. Vectorized meshgrid to get all absolute neighbor coordinates
            # absolute_coords shape: (N_t, M_t, 2)
            absolute_coords = type_points[:, np.newaxis, :] + rel_neighs[np.newaxis, :, :]
            
            # 3. Bulk Query (KDTree supports multidimensional input)
            # distances/indices shape: (N_t, M_t)
            distances, indices = self.kdtree.query(absolute_coords, k=1)

            # 4. Bulk Assignment
            if self.pbc:
                if np.any(distances > 1e-2 * self.toll_scale):
                    raise ValueError(f"Distance exceeds tolerance for type {t}")
                
                bigger_indices_arr[mask, :n_curr] = indices
                smaller_indices_arr[mask, :n_curr] = map_arr[indices]
            else:
                bigger_indices_arr[mask, :n_curr] = indices
                smaller_indices_arr[mask, :n_curr] = indices
                
            smaller_distances_arr[mask, :n_curr] = distances

        return bigger_indices_arr, smaller_indices_arr, smaller_distances_arr
    
    def get_neighbors_within_radius(self, query_points: np.ndarray, radius: float):
        assert self.kdtree is not None
        
        # 1. Scipy's fast radius search
        # Returns a list of lists: [[neighs_for_p0], [neighs_for_p1], ...]
        indices_list = self.kdtree.query_ball_point(query_points, r=radius)

        # 2. Vectorized Flattening
        # Calculate counts per query point to repeat indices correctly
        counts = np.array([len(nb) for nb in indices_list])
        if counts.sum() == 0:
            return np.array([], dtype=int), np.array([], dtype=int), np.array([]).reshape(0, 2)

        # query_indices: [0, 0, 0, 1, 1, ...]
        query_indices = np.repeat(np.arange(len(query_points)), counts)
        # tree_indices: flattened neighbor indices from the KDTree
        tree_indices = np.concatenate([nb for nb in indices_list if len(nb) > 0])

        # 3. Retrieve Coordinates
        if self.pbc:
            neighbor_coords = self.bigger_points[tree_indices]
            # REPLACEMENT FOR LAMBDA MAPPER: Direct Array Lookup
            # We pre-mapped these during generate_kdtree
            map_arr = np.array([self.mappings[i] for i in range(len(self.bigger_points))])
            lattice_indices = map_arr[tree_indices]
        else:
            neighbor_coords = self.points[tree_indices]
            lattice_indices = tree_indices

        return query_indices, lattice_indices, neighbor_coords

    # def plot_lattice(self, plot_connections: bool = True, colours: list = ["r", "g", "b", "c", "m", "y", "k"]) -> None:
    #     """
    #     Plots the lattice points and optionally the connections between them and the unit cell structure.

    #     Args:
    #         plot_connections (bool, optional): If True, plots the connections between nearest neighbors.
    #             Defaults to True.
    #         colours (list, optional): List of matplotlib colours to use for different point types. 

    #     Behavior:
    #         - Plots all lattice points grouped by type.
    #         - If `plot_connections` is True, it draws dashed red lines between nearest neighbors.

    #     Example:
    #     ```python
    #     lattice = MyLattice()
    #     lattice.generate_points()
    #     lattice.plot_lattice(plot_connections=True)
    #     ```

    #     Visualization Details:
    #         - Lattice points are plotted as small dots.
    #         - Nearest neighbor connections (if enabled) are shown as dashed red lines.
    #     """

    #     if len(colours) == 1: cols = {t[-1]:colours[0] for i, t, in enumerate(self.lattice_points)}
    #     else: cols = {t[-1]:colours[i] for i, t, in enumerate(self.lattice_points)}

    #     # plt.scatter(
    #     #     [self.points[:, 0]],
    #     #     [self.points[:, 1]],
    #     #     s=10, c=np.vectorize(cols.get)(self.point_types)
    #     # )

    #     if plot_connections:
    #         a = []
    #         for this, point_type in zip(self.points, self.point_types):
    #             for delta in self.neighbours[point_type]:
    #                 a.append([this[0], this[0] + delta[0]])
    #                 a.append((this[1], this[1] + delta[1]))
    #                 a.append("k--")
    #         # a = a[:30]
    #         # print(a)
    #         plt.plot(*a, alpha=0.1)

    #     # plt.title("Lattice Points")
    #     # plt.xlabel("X Coordinate")
    #     # plt.ylabel("Y Coordinate")
    #     plt.axis("equal")

    def plot_lattice(self, plot_connections: bool = True, colours: list = ["r", "g", "b", "c", "m", "y", "k"]) -> None:
        # Create a color map: index -> color
        # Since self.point_type_names gives ["A", "B", ...], we map 0 -> 'r', 1 -> 'g', etc.
        if len(colours) == 1:
            color_map = {i: colours[0] for i in range(len(self.point_type_names))}
        else:
            color_map = {i: colours[i % len(colours)] for i in range(len(self.point_type_names))}

        # 1. Plot the points themselves
        # We use a vectorized approach for the colors based on the integer point_types
        point_colors = [color_map[t] for t in self.point_types]
        plt.scatter(self.points[:, 0], self.points[:, 1], s=10, c=point_colors, zorder=3)

        # 2. Plot connections (hopping lines)
        if plot_connections:
            lines = []
            # self.neighbours is a dict: {"A": [[x,y], ...]}
            # self.point_type_names is a list: ["A", "B"]
            for i, (this_pos, type_idx) in enumerate(zip(self.points, self.point_types)):
                name = self.point_type_names[type_idx]
                for delta in self.neighbours[name]:
                    lines.append([this_pos[0], this_pos[0] + delta[0]]) # x-coords
                    lines.append([this_pos[1], this_pos[1] + delta[1]]) # y-coords
                    lines.append("k--")
            
            if lines:
                plt.plot(*lines, alpha=0.1, zorder=2)

        plt.axis("equal")


    def __repr__(self):
        return (
            f"""Layer(
    lv1 = {self.lv1},
    lv2 = {self.lv2},
    lattice_points = {self.lattice_points},
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
    def pbc(self):
        return self._rust_lattice.pbc
    @pbc.setter
    def pbc(self, value):
        raise LatticeAlreadyFinalisedError("pbc", self.__class__.__name__)

    @property
    def toll_scale(self):
        return self._rust_lattice.toll_scale
    @toll_scale.setter
    def toll_scale(self, value):
        raise LatticeAlreadyFinalisedError("toll_scale", self.__class__.__name__)

    @property
    def study_proximity(self):
        return self._rust_lattice.study_proximity
    @study_proximity.setter
    def study_proximity(self, value):
        raise LatticeAlreadyFinalisedError("study_proximity", self.__class__.__name__)

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
    def point_type_names(self):
        return self._rust_lattice.point_type_names
    @point_type_names.setter
    def point_type_names(self, value):
        raise LatticeAlreadyFinalisedError("point_type_names", self.__class__.__name__)

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
    def bigger_points(self):
        return self._rust_lattice.bigger_points
    @bigger_points.setter
    def bigger_points(self, value):
        raise LatticeAlreadyFinalisedError("bigger_points", self.__class__.__name__)

    @property
    def mappings(self):
        return self._rust_lattice.mappings
    @mappings.setter
    def mappings(self, value):
        raise LatticeAlreadyFinalisedError("mappings", self.__class__.__name__)

    @property
    def neighbours(self):
        return self._rust_lattice.neighbours
    @neighbours.setter
    def neighbours(self, value):
        raise LatticeAlreadyFinalisedError("neighbours", self.__class__.__name__)


# ===============================================
# ============= some example layers =============
# ===============================================

class SquareLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        # local definitions, passed to parent
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0, 1])  # Lattice vector in the y-direction
        lattice_points = [
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
        super().__init__(lv1, lv2, lattice_points, neighbours, pbc, study_proximity)



class TriangularLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3) / 2])  # Lattice vector at 60 degrees
        lattice_points = [
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
        super().__init__(lv1, lv2, lattice_points, neighbours, pbc, study_proximity)



class Rhombus60Layer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        angle = 60  # hardcoded angle... make a copy of the whole class for different angles
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        cos_angle = np.cos(np.radians(angle))
        sin_angle = np.sin(np.radians(angle))
        lv2 = np.array([cos_angle, sin_angle])  # Lattice vector at specified angle
        lattice_points = [
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
        super().__init__(lv1, lv2, lattice_points, neighbours, pbc, study_proximity)



class KagomeLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0])  # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3)/2])  # Lattice vector at 60 degrees

        lattice_points = [
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
        super().__init__(lv1, lv2, lattice_points, neighbours, pbc, study_proximity)


class HexagonalLayer(Layer):
    def __init__(self, pbc=False, study_proximity: int=1) -> None:
        lv1 = np.array([1, 0]) # Lattice vector in the x-direction
        lv2 = np.array([0.5, np.sqrt(3) / 2])

        lattice_points = [
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
        super().__init__(lv1, lv2, lattice_points, neighbours, pbc, study_proximity)

