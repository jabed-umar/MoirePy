import numpy as np
from moirepy.presets import SquareLayer
from moirepy.stack import MoireStack

def test_phase():
    l1 = SquareLayer()
    stack = MoireStack(l1)
    stack.generate_points(1, 1)
    
    def t_intra(coo_i, coo_j, R, labels_i, labels_j, stack):
        return np.ones(len(coo_i))
        
    stack.generate_hamiltonian(t_matrix={}, t_self={0: t_intra})
    phase = stack.get_phase([0.1, 0.2])
    assert len(phase) == len(stack.H.data)

