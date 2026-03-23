import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from utils import OUT_JSON, OUT_EXP1, OUT_EXP2

def plot_results(x: np.ndarray, y: np.ndarray, out_path):
    plt.rcParams.update({
        "font.size": 14, "axes.labelsize": 18, "xtick.labelsize": 14,
        "ytick.labelsize": 14, "lines.linewidth": 2.0, "lines.markersize": 6,
    })

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
    
    # Save as PDF as requested
    final_path = out_path.with_suffix(".pdf")
    plt.savefig(final_path, bbox_inches="tight")
    
    # # save webp for faster loading in websites
    # final_path = out_path.with_suffix(".webp")
    # plt.savefig(final_path, bbox_inches="tight", format="webp", pil_kwargs={"quality": 60, "lossless": False, "method": 6, "optimize": True})

    plt.close()
    print(f"Plot saved to {final_path}")

def main():
    if not OUT_JSON.exists():
        print(f"Error: {OUT_JSON} not found. Run benchmark.py first.")
        return

    with open(OUT_JSON, "r") as f:
        data = json.load(f)

    results = data["results"]
    x = np.array([r["cells_times_2"] for r in results])
    y1 = np.array([r["experiment_1_ms"] for r in results])
    y2 = np.array([r["experiment_2_ms"] for r in results])

    plot_results(x, y1, OUT_EXP1)
    plot_results(x, y2, OUT_EXP2)

if __name__ == "__main__":
    main()