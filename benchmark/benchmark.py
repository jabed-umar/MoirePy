import argparse
import csv
import json
from pathlib import Path
import platform
import subprocess
import tarfile
from time import time

import numpy as np
from tqdm import trange

from moirepy import (
    BilayerMoireLattice,
    HexagonalLayer,
    KagomeLayer,
    Rhombus60Layer,
    SquareLayer,
    TriangularLayer,
)


LATTICE_TYPES = {
    "hexagonal": HexagonalLayer,
    "square": SquareLayer,
    "triangular": TriangularLayer,
    "rhombus60": Rhombus60Layer,
    "kagome": KagomeLayer,
}

LATTICE_TYPE_KEY = "hexagonal"
N1 = 1
N2 = 1
INTER_LAYER_RADIUS = 1.0
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_N = 25
DEFAULT_REPEAT = 100
DEFAULT_DATA_PATH = SCRIPT_DIR / "commensurate_angles.csv"
DEFAULT_OUT_DIR = SCRIPT_DIR
OUTPUT_JSON_NAME = "benchmark_results.json"
HOPPING = {
    "tll": 1.0,
    "tuu": 1.0,
    "tlu": 1.0,
    "tul": 1.0,
    "tlself": 1.0,
    "tuself": 1.0,
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MoirePy benchmark and save JSON results.")

    parser.add_argument("--data-path", type=Path, default=DEFAULT_DATA_PATH)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)

    parser.add_argument("--n", type=int, default=DEFAULT_N)
    parser.add_argument("--repeat", type=int, default=DEFAULT_REPEAT)
    parser.add_argument("--cpu_name", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")

    return parser.parse_args()


def first_non_empty_line(text: str) -> str | None:
    for line in text.splitlines():
        candidate = line.strip()
        if not candidate:
            continue
        if candidate.lower() == "name":
            continue
        return candidate
    return None


def run_command_output(command: list[str]) -> str | None:
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=2,
            check=False,
        )
    except Exception:
        return None

    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None


def get_cpu_name() -> str | None:
    system = platform.system().lower()

    if system == "linux":
        cpuinfo = Path("/proc/cpuinfo")
        if cpuinfo.exists():
            try:
                for line in cpuinfo.read_text(encoding="utf-8", errors="ignore").splitlines():
                    if line.lower().startswith("model name"):
                        _, value = line.split(":", 1)
                        cpu_name = value.strip()
                        if cpu_name:
                            return cpu_name
            except Exception:
                pass

        output = run_command_output(["lscpu"])
        if output:
            for line in output.splitlines():
                if "model name:" in line.lower():
                    _, value = line.split(":", 1)
                    cpu_name = value.strip()
                    if cpu_name:
                        return cpu_name

    elif system == "darwin":
        output = run_command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
        if output:
            cpu_name = first_non_empty_line(output)
            if cpu_name:
                return cpu_name

    elif system == "windows":
        output = run_command_output(
            [
                "powershell",
                "-NoProfile",
                "-Command",
                "(Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)",
            ]
        )
        if output:
            cpu_name = first_non_empty_line(output)
            if cpu_name:
                return cpu_name

        output = run_command_output(["wmic", "cpu", "get", "name"])
        if output:
            cpu_name = first_non_empty_line(output)
            if cpu_name:
                return cpu_name

    fallback = platform.processor().strip() if platform.processor() else ""
    if fallback:
        return fallback

    fallback = platform.uname().processor.strip() if platform.uname().processor else ""
    if fallback:
        return fallback

    return None


def extract_tar_gz(tar_path: Path, dest_dir: Path) -> None:
    with tarfile.open(tar_path, "r:gz") as tar:
        try:
            tar.extractall(path=dest_dir, filter="data")
        except TypeError:
            tar.extractall(path=dest_dir)


def ensure_csv_exists(csv_path: Path) -> Path:
    if csv_path.exists():
        return csv_path

    tar_path = csv_path.with_name(f"{csv_path.stem}.tar.gz")
    if not tar_path.exists():
        raise FileNotFoundError(f"Data file not found: {csv_path}")

    print(f"CSV not found. Trying to extract from {tar_path} ...")
    extract_tar_gz(tar_path, csv_path.parent)

    if csv_path.exists():
        return csv_path

    for candidate in csv_path.parent.rglob(csv_path.name):
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        f"Data file not found after extracting tar.gz: {csv_path} (from {tar_path})"
    )


def load_rows(csv_path: Path) -> list[dict[str, float | int]]:
    csv_path = ensure_csv_exists(csv_path)

    rows: list[dict[str, float | int]] = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "idx": int(row.get("", "0") or 0),
                    "angle_deg": float(row["angle (deg)"]),
                    "angle_rad": float(row["angle (rad)"]),
                    "ll1": int(row["ll1"]),
                    "ll2": int(row["ll2"]),
                    "ul1": int(row["ul1"]),
                    "ul2": int(row["ul2"]),
                    "cells": int(round(float(row["cells"]))),
                }
            )
    return rows


def print_run_configuration(args: argparse.Namespace, cpu_name: str, cpu_source: str) -> None:
    print("Run configuration")
    print("-----------------")

    if str(args.data_path) == str(DEFAULT_DATA_PATH):
        print(f"Data Path: {args.data_path} (default; change with --data-path)")
    else:
        print(f"Data Path: {args.data_path}")

    if str(args.out_dir) == str(DEFAULT_OUT_DIR):
        print(f"Output Directory: {args.out_dir} (default; change with --out-dir)")
    else:
        print(f"Output Directory: {args.out_dir}")

    print(f"Number of points (n): {args.n}")
    print(f"Repeats: {args.repeat}")

    if cpu_source == "auto-detected":
        print(f"cpu_name: {cpu_name} (auto-detected, override with --cpu_name)")
    elif cpu_source == "manual":
        print(f"cpu_name: {cpu_name} (manual via --cpu_name)")
    else:
        print(f'cpu_name: {cpu_name} (auto-detect failed, use --cpu_name "<your cpu name>")')

    print()


def select_equally_spaced_rows(rows: list[dict[str, float | int]], n: int) -> list[dict[str, float | int]]:
    if n <= 0:
        raise ValueError("n must be positive")
    if len(rows) < n:
        raise ValueError(f"Not enough rows in data ({len(rows)}) for n={n}")

    sorted_rows = sorted(rows, key=lambda r: int(r["cells"]))
    targets = np.linspace(int(sorted_rows[0]["cells"]), int(sorted_rows[-1]["cells"]), n)

    remaining = sorted_rows.copy()
    selected: list[dict[str, float | int]] = []

    for target in targets:
        i = min(range(len(remaining)), key=lambda j: abs(int(remaining[j]["cells"]) - target))
        selected.append(remaining.pop(i))

    return sorted(selected, key=lambda r: int(r["cells"]))


def build_lattice(row: dict[str, float | int], lattice_class, n1: int, n2: int) -> BilayerMoireLattice:
    return BilayerMoireLattice(
        lattice_class,
        int(row["ll1"]),
        int(row["ll2"]),
        int(row["ul1"]),
        int(row["ul2"]),
        n1,
        n2,
        verbose=False,
    )


def run_experiment_1(
    row: dict[str, float | int],
    repeat: int,
    lattice_class,
    n1: int,
    n2: int,
    inter_layer_radius: float,
    hopping: dict[str, float],
    quiet: bool,
) -> float:
    start = time()
    for _ in trange(repeat, desc=f"Exp1 cells={row['cells']}", leave=False, disable=quiet):
        lattice = build_lattice(row, lattice_class, n1, n2)
        lattice.generate_connections(inter_layer_radius=inter_layer_radius)
        lattice.generate_hamiltonian(**hopping)
    return (time() - start) * 1000.0 / repeat


def run_experiment_2(
    row: dict[str, float | int],
    repeat: int,
    lattice_class,
    n1: int,
    n2: int,
    inter_layer_radius: float,
    hopping: dict[str, float],
    quiet: bool,
) -> float:
    lattice = build_lattice(row, lattice_class, n1, n2)
    lattice.generate_connections(inter_layer_radius=inter_layer_radius)

    start = time()
    for _ in trange(repeat, desc=f"Exp2 cells={row['cells']}", leave=False, disable=quiet):
        lattice.generate_hamiltonian(**hopping)
    return (time() - start) * 1000.0 / repeat


def main() -> None:
    args = parse_args()
    lattice_class = LATTICE_TYPES[LATTICE_TYPE_KEY]
    out_json_path = args.out_dir / OUTPUT_JSON_NAME
    cpu_name = args.cpu_name.strip() if args.cpu_name else None
    cpu_source = "manual"

    if cpu_name is None:
        cpu_name = get_cpu_name()
        if cpu_name is None:
            cpu_name = "unknown"
            cpu_source = "failed"
            print('Tried to fetch CPU name, but failed. Run again with --cpu_name "<your cpu name>".')
        else:
            cpu_source = "auto-detected"

    print_run_configuration(args, cpu_name, cpu_source)

    rows = load_rows(args.data_path)
    selected_rows = select_equally_spaced_rows(rows, args.n)

    if not args.quiet:
        print(f"Loaded rows: {len(rows)}")
        print(f"Selected rows: {len(selected_rows)}")
        print("cells | ll1 ll2 ul1 ul2")
        for row in selected_rows:
            print(f"{row['cells']} | {row['ll1']} {row['ll2']} {row['ul1']} {row['ul2']}")
        print()

    x_values = np.array([int(row["cells"]) * 2 for row in selected_rows], dtype=float)

    exp1_times = np.array(
        [
            run_experiment_1(
                row,
                args.repeat,
                lattice_class,
                N1,
                N2,
                INTER_LAYER_RADIUS,
                HOPPING,
                args.quiet,
            )
            for row in selected_rows
        ],
        dtype=float,
    )

    exp2_times = np.array(
        [
            run_experiment_2(
                row,
                args.repeat,
                lattice_class,
                N1,
                N2,
                INTER_LAYER_RADIUS,
                HOPPING,
                args.quiet,
            )
            for row in selected_rows
        ],
        dtype=float,
    )

    payload = {
        "config": {
            "n": args.n,
            "repeat": args.repeat,
            "cpu_name": cpu_name,
            "platform_system": platform.system(),
            "platform_machine": platform.machine(),
            "x_axis": "cells * 2",
            "y_axis": "time taken (ms), averaged over repeat runs",
        },
        "results": [
            {
                "cells_times_2": int(x),
                "experiment_1_ms": float(t1),
                "experiment_2_ms": float(t2),
            }
            for row, x, t1, t2 in zip(selected_rows, x_values, exp1_times, exp2_times)
        ],
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Benchmark complete. Results saved to {out_json_path}")


if __name__ == "__main__":
    main()
