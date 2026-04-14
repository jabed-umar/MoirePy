import pytest
import os

# --- The Bouncer (Whitelist) ---
ALLOWED_TESTS = {
    "quickstart_tbg.py",
    "nori.py",
    
    # "custom_hoppings.py",
    # "k_space_hamiltonian.py",
    # "multiple_orbitals.py",
    # "non_hermitian.py",
    # "obc_vs_pbc.py",
    
    # "koshino.py",  # too slow for CI, but will keep the tests writen, when I have time later 
    # "prepare_layers.py",  # nothing to test
}

def pytest_ignore_collect(collection_path, config):
    """Tell pytest which files to ignore."""
    # Ensure we use absolute paths to compare reliably
    abs_path = os.path.abspath(str(collection_path))
    HERE = os.path.dirname(os.path.abspath(__file__))
    
    # We only care about files inside the 'tests/' directory
    if abs_path.startswith(HERE):
        name = os.path.basename(abs_path)
        
        # 1. Never ignore conftest.py or the data folder itself
        if name in ("conftest.py", "data", "__init__.py"):
            return False
            
        # 2. If it's a file, check if it's on the whitelist
        if os.path.isfile(abs_path):
            if name not in ALLOWED_TESTS:
                return True # Kick it out
                
    return False

# --- CLI Options ---
def pytest_addoption(parser):
    parser.addoption(
        "--record", action="store_true", default=False, help="Record ground truth data"
    )

@pytest.fixture
def record_mode(request):
    return request.config.getoption("--record")
