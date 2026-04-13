import pytest
import numpy as np
import os
from testbook import testbook

# Absolute-ish paths relative to this file
HERE = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.abspath(os.path.join(HERE, "..", "docs", "getting_started", "quickstart_tbg.ipynb"))
DATA_DIR = os.path.join(HERE, "data", "quickstart_tbg")
EVALS_FILE = os.path.join(DATA_DIR, "evals.bin")

@testbook(NB_PATH, execute=True)
def test_quickstart_tbg_eigenvalues(tb, record_mode):
    tb.inject("""
import numpy as np
_test_dense_ham = k_hamiltonian.toarray()
_test_evals = np.linalg.eigvalsh(_test_dense_ham)
_test_evals_list = _test_evals.tolist()
    """)

    current_evals = np.array(tb.value("_test_evals_list"), dtype=np.float16)

    if record_mode:
        os.makedirs(DATA_DIR, exist_ok=True)
        current_evals.tofile(EVALS_FILE)
        pytest.skip(f"Recorded to {EVALS_FILE}")
    else:
        if not os.path.exists(EVALS_FILE):
            pytest.fail("Reference data not found. Run with --record.")
            
        ref_evals = np.fromfile(EVALS_FILE, dtype=np.float16)
        
        if len(current_evals) != len(ref_evals):
            pytest.fail(f"Shape mismatch: {len(current_evals)} vs {len(ref_evals)}")
            
        assert np.allclose(current_evals, ref_evals)
