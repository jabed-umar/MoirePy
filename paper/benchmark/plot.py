import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate benchmark plots from a benchmark JSON file.")
    parser.add_argument("--input-json", type=Path, required=True, help="Path to benchmark JSON file.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory where PDF plots will be written.",
    )
    return parser.parse_args()


def plot_results(x: np.ndarray, y: np.ndarray, out_path: Path) -> None:
    plt.rcParams.update(
        {
            "font.size": 14,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "lines.linewidth": 2.0,
            "lines.markersize": 6,
        }
    )

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(x, y, marker="o")
    ax.set_xlabel("Number of lattice points")
    ax.set_ylabel("Time taken (ms)")

    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((3, 3))
    ax.xaxis.set_major_formatter(formatter)

    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to {out_path}")


def main() -> None:
    args = parse_args()

    if not args.input_json.exists():
        raise FileNotFoundError(f"Input JSON not found: {args.input_json}")

    data = json.loads(args.input_json.read_text(encoding="utf-8"))
    results = data.get("results", [])

    if not results:
        raise ValueError(f"No results found in JSON: {args.input_json}")

    x = np.array([r["cells_times_2"] for r in results], dtype=float)
    y1 = np.array([r["experiment_1_ms"] for r in results], dtype=float)
    y2 = np.array([r["experiment_2_ms"] for r in results], dtype=float)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    out_exp1 = args.output_dir / "experiment1_total_pipeline.pdf"
    out_exp2 = args.output_dir / "experiment2_hamiltonian_only.pdf"

    plot_results(x, y1, out_exp1)
    plot_results(x, y2, out_exp2)


if __name__ == "__main__":
    main()
