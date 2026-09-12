import numpy as np
from moirepy.presets import SquareLayer

def test_neighbors():
    layer = SquareLayer()
    layer.generate_points(np.array([2.0, 0.0, 0.0]), np.array([0.0, 2.0, 0.0]))
    
    i, j, li, lj = layer.first_nearest_neighbours()
    assert len(i) == 4 * 4 # 4 points, each has 4 neighbors
    
