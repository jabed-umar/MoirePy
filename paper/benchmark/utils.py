from pathlib import Path
import csv
import numpy as np
from moirepy import HexagonalLayer

# Benchmark controls
N = 25
REPEAT = 100
LATTICE_TYPE = HexagonalLayer

# Paths
folder = ""
# folder = "pi5/"
# folder = "i3_gen4/"
DATA_PATH = Path(__file__).parent / f"{folder}commensurate_angles.csv"
OUT_EXP1 = Path(__file__).parent / f"{folder}experiment1_total_pipeline"
OUT_EXP2 = Path(__file__).parent / f"{folder}experiment2_hamiltonian_only"
OUT_JSON = Path(__file__).parent / f"{folder}benchmark_results.json"

def load_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    if not path.exists():
        raise FileNotFoundError(f"Data file not found at {path}")
        
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "idx": int(r.get("", "0") or 0),
                "angle_deg": float(r["angle (deg)"]),
                "angle_rad": float(r["angle (rad)"]),
                "ll1": int(r["ll1"]),
                "ll2": int(r["ll2"]),
                "ul1": int(r["ul1"]),
                "ul2": int(r["ul2"]),
                "cells": int(round(float(r["cells"]), 0)),
            })
    return rows

def select_equally_spaced_by_cells(rows: list[dict[str, float | int]], n: int) -> list[dict[str, float | int]]:
    if n <= 0: raise ValueError("N must be positive")
    if len(rows) < n: raise ValueError(f"Not enough rows in data ({len(rows)}) for N={n}")

    sorted_rows = sorted(rows, key=lambda r: int(r["cells"]))
    min_cells = int(sorted_rows[0]["cells"])
    max_cells = int(sorted_rows[-1]["cells"])

    targets = np.linspace(min_cells, max_cells, n)
    remaining = sorted_rows.copy()
    picked: list[dict[str, float | int]] = []

    for t in targets:
        i = min(range(len(remaining)), key=lambda j: abs(int(remaining[j]["cells"]) - t))
        picked.append(remaining.pop(i))

    return sorted(picked, key=lambda r: int(r["cells"]))

def print_head(rows: list[dict[str, float | int]]) -> None:
    print("Head (schema check):")
    print("cells | angle (deg) | ll1 | ll2 | ul1 | ul2 |")
    for r in rows:
        print(f"{r['cells']} | {r['angle_deg']:.4f} | {r['ll1']} | {r['ll2']} | {r['ul1']} | {r['ul2']}")
    print()
