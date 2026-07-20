import numpy as np
from moirepy.layers import Layer

def test_generate_points_square():
    layer = Layer([1.0, 0.0], [0.0, 1.0], [(0,0,"A")], {})
    layer.generate_points([2.0, 0.0], [0.0, 2.0])
    
    assert len(layer._points) == 4
    assert len(layer._bigger_points) == 36
    assert len(layer._bigger_to_primary) == 36

def test_generate_points_rotated():
    layer = Layer([1.0, 0.0], [0.0, 1.0], [(0,0,"A")], {})
    layer.perform_rotation_translation(90.0)
    layer.generate_points([2.0, 0.0], [0.0, 2.0])
    
    assert len(layer._points) == 4
    assert len(layer._bigger_points) == 36
