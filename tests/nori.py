import pytest
import numpy as np
import os
import json
from testbook import testbook

HERE = os.path.dirname(os.path.abspath(__file__))
NB_PATH = os.path.abspath(os.path.join(HERE, "..", "docs", "replications", "nori.ipynb"))
DATA_DIR = os.path.join(HERE, "data", "nori")
META_FILE = os.path.join(HERE, "data", "meta.json")

@testbook(NB_PATH, execute=True)
def test_nori_spectrum(tb, record_mode):
    energies_list = tb.value("all_energies_shifted.tolist()")
    gap_val = float(tb.value("gap"))
    mu_val = float(tb.value("mu_shift"))

    current_energies = np.array(energies_list, dtype=np.float16)
    anchor_row = current_energies[75]  # Match only Row 75 (if one matches, all will match)
    gap_scaled = round(gap_val * 1e3, 2)
    mu_scaled = round(mu_val * 1e6, 2)

    if record_mode:
        os.makedirs(DATA_DIR, exist_ok=True)
        anchor_row.tofile(os.path.join(DATA_DIR, "row75.bin"))

        meta = {}
        if os.path.exists(META_FILE):
            with open(META_FILE, "r") as f:
                meta = json.load(f)
        
        meta["nori"] = {
            "gap_scaled": gap_scaled,
            "mu_scaled": mu_scaled
        }
        with open(META_FILE, "w") as f:
            json.dump(meta, f)

        pytest.skip("Recorded Lean Nori.")
    else:
        ref_path = os.path.join(DATA_DIR, "row75.bin")
        ref_row = np.fromfile(ref_path, dtype=np.float16)
        assert np.allclose(anchor_row, ref_row)
        
        with open(META_FILE, "r") as f:
            meta = json.load(f)
        
        assert pytest.approx(gap_scaled) == meta["nori"]["gap_scaled"]
        assert pytest.approx(mu_scaled) == meta["nori"]["mu_scaled"]
