import numpy as np
import pytest
from moirepy.layers import Layer
from moirepy.utils import LatticeAlreadyFinalisedError

def test_layer_init():
    lv1 = [1.0, 0.0]
    lv2 = [0.5, 0.866]
    basis = [(0, 0, "A"), (0.5, 0.288, 0.1, "B")]
    
    layer = Layer(lv1, lv2, basis, {})
    
    assert layer.lv1.shape == (3,)
    assert layer.lv2.shape == (3,)
    assert np.allclose(layer.basis_points[1], [0.5, 0.288, 0.1])
    assert layer.orbitals == 1

def test_layer_strain():
    lv1 = [1.0, 0.0]
    lv2 = [0.0, 1.0]
    layer = Layer(lv1, lv2, [(0,0,"A")], {})
    
    # 10% stretch in x
    strain = [[0.1, 0.0], [0.0, 0.0]]
    layer.apply_strain(strain)
    
    assert np.allclose(layer.lv1, [1.1, 0.0, 0.0])
    assert np.allclose(layer.lv2, [0.0, 1.0, 0.0])

def test_layer_strain_realignment():
    lv1 = [1.0, 0.0]
    lv2 = [0.0, 1.0]
    layer = Layer(lv1, lv2, [(0,1,"A")], {})
    
    # Shear strain that rotates lv1
    strain = [[0.0, 0.0], [0.1, 0.0]]
    layer.apply_strain(strain)
    
    # lv1 should still be along x-axis!
    assert np.isclose(layer.lv1[1], 0.0)
    assert np.allclose(layer.lv2, [0.09950372, 0.99503719, 0.0], atol=1e-3)

def test_layer_validation():
    with pytest.raises(ValueError):
        Layer([1.0, 0.1], [0.0, 1.0], [], {}) # lv1 not on x-axis
    
    with pytest.raises(ValueError):
        Layer([1.0, 0.0], [1.0, -0.1], [], {}) # lv2 y < 0
