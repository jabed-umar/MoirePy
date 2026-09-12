import numpy as np
from moirepy.presets import SquareLayer
from moirepy.stack import MoireStack

def test_connections():
    l1 = SquareLayer()
    l2 = SquareLayer()
    stack = MoireStack(l1)
    stack.add_layer(l2, 1, 0, 1, 0, translate=(0.5, 0.5, 3.0))
    stack.generate_points(1, 1)
    
    radii = np.array([[0, 4.0], [4.0, 0]])
    stack.generate_connections(radii)
    
    assert (0, 1) in stack.connections
    assert (1, 0) in stack.connections
