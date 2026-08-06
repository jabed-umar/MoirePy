import numpy as np
import pytest
from moirepy.presets import HexagonalLayer, SquareLayer
from moirepy.stack import MoireStack

def test_stack_basic():
    l1 = HexagonalLayer()
    stack = MoireStack(l1)
    
    assert stack.state == "BUILDING"
    assert stack.mlv1 is None
    
    l2 = HexagonalLayer()
    # 21.78 deg TBG (ll=(1,2), ul=(2,1))
    stack.add_layer(l2, 1, 2, 2, 1)
    
    assert stack.mlv1 is not None
    assert stack.mlv2 is not None
    assert len(stack.layers) == 2
    
def test_stack_commensurate_fail():
    l1 = HexagonalLayer()
    stack = MoireStack(l1)
    
    l2 = HexagonalLayer()
    # Apply strain to break commensurability
    l2.apply_strain([[0.1, 0.0], [0.0, 0.0]])
    
    with pytest.raises(ValueError, match="Commensurate length mismatch"):
        stack.add_layer(l2, 1, 2, 2, 1)

def test_stack_hetero_commensurate():
    l1 = HexagonalLayer()
    stack = MoireStack(l1)
    
    # Fake hetero layer, starting with different lattice constant
    l2 = HexagonalLayer()
    l2.lv1 *= 1.05
    l2.lv2 *= 1.05
    
    # We want to form same supercell, so we must strain l2 back!
    # If we don't, it will fail
    with pytest.raises(ValueError):
        stack.add_layer(l2, 1, 2, 2, 1)
        
    # We apply strain to make it perfectly match
    l3 = HexagonalLayer()
    l3.lv1 *= 1.05
    l3.lv2 *= 1.05
    
    # We know we need 1/1.05 = 0.95238...
    s = (1.0 / 1.05) - 1.0
    l3.apply_strain([[s, 0.0], [0.0, s]])
    
    stack.add_layer(l3, 1, 2, 2, 1)
    assert len(stack.layers) == 2

def test_stack_generate_points():
    l1 = HexagonalLayer()
    stack = MoireStack(l1)
    
    l2 = HexagonalLayer()
    stack.add_layer(l2, 1, 2, 2, 1)
    
    assert stack.state == "BUILDING"
    stack.generate_points()
    assert stack.state == "GENERATED"
    
    assert len(l1._points) > 0
    assert len(l2._points) > 0
    
    l3 = HexagonalLayer()
    with pytest.raises(Exception):
        stack.add_layer(l3, 1, 2, 2, 1)

def test_stack_repeat():
    l1 = HexagonalLayer()
    stack = MoireStack(l1)
    
    stack.generate_points()
    n_pts_single = len(l1._points)
    
    stack.repeat(2, 2)
    assert len(l1._points) == n_pts_single * 4
    
    stack.repeat(1, 2)
    assert len(l1._points) == n_pts_single * 8
