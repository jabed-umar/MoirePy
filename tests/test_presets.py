import numpy as np
from moirepy.presets import SquareLayer, HexagonalLayer

def test_square_preset():
    layer = SquareLayer()
    assert np.allclose(layer.lv1, [1.0, 0.0, 0.0])
    assert layer.basis_points.shape == (1, 3)

def test_hexagonal_preset():
    layer = HexagonalLayer()
    assert np.allclose(layer.lv2, [0.5, np.sqrt(3.0)/2.0, 0.0])
    assert layer.basis_points.shape == (2, 3)
