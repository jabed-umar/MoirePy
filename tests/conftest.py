import pytest
import os

# --- The Bouncer (Whitelist) ---
ALLOWED_TESTS = {
    "test_quickstart_tbg.py",
}

def pytest_ignore_collect(collection_path, config):
    """Tell pytest which files to ignore."""
    rel_path = os.path.relpath(str(collection_path), os.getcwd())
    
    # We only care about files inside the 'tests/' directory
    if rel_path.startswith("tests"):
        name = os.path.basename(rel_path)
        
        # 1. Never ignore conftest.py or the data folder itself
        if name in ("conftest.py", "data", "__init__.py"):
            return False

        # 2. If it's a file, check if it's on the whitelist
        if os.path.isfile(collection_path):
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
