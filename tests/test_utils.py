import numpy as np
import pytest
from moirepy.utils import get_rotation_matrix, are_coeffs_integers, LatticeAlreadyFinalisedError, StackStateError

def test_get_rotation_matrix():
    rot = get_rotation_matrix(np.pi / 2)
    # Rotating (1, 0, 5) by 90 degrees around Z should give (0, 1, 5)
    vec = np.array([1.0, 0.0, 5.0])
    res = rot @ vec
    np.testing.assert_allclose(res, [0.0, 1.0, 5.0], atol=1e-10)

def test_are_coeffs_integers():
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([0.0, 1.0, 0.0])
    
    # Integer combination
    v3 = 3 * v1 - 2 * v2
    assert are_coeffs_integers(v1, v2, v3) is True
    
    # Non-integer combination
    v4 = 1.5 * v1 + v2
    assert are_coeffs_integers(v1, v2, v4) is False

def test_exceptions():
    with pytest.raises(LatticeAlreadyFinalisedError):
        raise LatticeAlreadyFinalisedError("test")
        
    with pytest.raises(StackStateError):
        raise StackStateError("test")
