import json
from time import time
import numpy as np
from tqdm import trange
from moirepy import BilayerMoireLattice

from utils import (
    N, REPEAT, LATTICE_TYPE, DATA_PATH, OUT_JSON, OUT_EXP1, OUT_EXP2,
    load_rows, select_equally_spaced_by_cells, print_head
)

def run_experiment_1(row: dict[str, float | int], repeat: int) -> float:
    ll1, ll2 = int(row["ll1"]), int(row["ll2"])
    ul1, ul2 = int(row["ul1"]), int(row["ul2"])
    t0 = time()
    for _ in trange(repeat, desc=f"Exp1 cells={row['cells']}", leave=False):
        lattice = BilayerMoireLattice(LATTICE_TYPE, ll1, ll2, ul1, ul2, 1, 1, verbose=False)
        lattice.generate_connections(inter_layer_radius=1)
        lattice.generate_hamiltonian(tll=1, tuu=1, tlu=1, tul=1, tuself=1, tlself=1)
    return (time() - t0) * 1000.0 / repeat

def run_experiment_2(row: dict[str, float | int], repeat: int) -> float:
    ll1, ll2 = int(row["ll1"]), int(row["ll2"])
    ul1, ul2 = int(row["ul1"]), int(row["ul2"])
    lattice = BilayerMoireLattice(LATTICE_TYPE, ll1, ll2, ul1, ul2, 1, 1, verbose=False)
    lattice.generate_connections(inter_layer_radius=1)
    t0 = time()
    for _ in trange(repeat, desc=f"Exp2 cells={row['cells']}", leave=False):
        lattice.generate_hamiltonian(tll=1, tuu=1, tlu=1, tul=1, tuself=1, tlself=1)
    return (time() - t0) * 1000.0 / repeat

def save_results_json(sampled, x_vals, y_exp1, y_exp2):
    payload = {
        "config": {"N": N, "REPEAT": REPEAT, "LATTICE_TYPE": LATTICE_TYPE.__name__},
        "results": [
            {
                "cells_times_2": int(x),
                "experiment_1_ms": float(t1),
                "experiment_2_ms": float(t2),
            }
            for x, t1, t2 in zip(x_vals, y_exp1, y_exp2)
        ],
        "artifacts": {"plot_1": str(OUT_EXP1), "plot_2": str(OUT_EXP2)}
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))

def main():
    rows = load_rows(DATA_PATH)
    sampled = select_equally_spaced_by_cells(rows, N)
    print_head(sampled)
    
    x_vals = np.array([int(r["cells"]) * 2 for r in sampled], dtype=float)
    y_exp1 = np.array([run_experiment_1(r, REPEAT) for r in sampled], dtype=float)
    y_exp2 = np.array([run_experiment_2(r, REPEAT) for r in sampled], dtype=float)

    save_results_json(sampled, x_vals, y_exp1, y_exp2)
    print(f"Benchmark complete. Results saved to {OUT_JSON}")

if __name__ == "__main__":
    main()