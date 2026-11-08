"""Plot saved or recomputed convergence-order tests; do not modify errors."""

from pathlib import Path
import argparse
import csv
import numpy as np
from scipy.io import loadmat
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--data-dir", type=Path, default=HERE / "data" / "order")
parser.add_argument("--output-dir", type=Path, default=HERE / "figures" / "order")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "font.size": 14,
        "axes.labelsize": 19,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
    }
)
fig, axes = plt.subplots(1, 2, figsize=(11.92, 4.8), layout="constrained")
rows = []
for ax, method, expected, color in zip(
    axes, ["CP", "N"], [1, 2], ["#0072BD", "#D95319"]
):
    r = loadmat(args.data_dir / f"d100_{method}.mat", simplify_cells=True)["r"]
    e = np.asarray(r["error"]).ravel()
    x, y = e[:-1], e[1:]
    assert np.all(e > 0) and r["converged"]
    # Predeclared local window excludes the initial transient and roundoff floor.
    local = (x < 1e-2) & (y > 1e-13) & (y < x)
    assert local.sum() >= 3
    slope, intercept = np.polyfit(np.log10(x[local]), np.log10(y[local]), 1)
    coefficient = np.median(y[local] / x[local] ** expected)
    indices = np.arange(0, len(x), 50 if method == "CP" else 1)
    ax.loglog(x[indices], y[indices], "o", color=color, markersize=4, label="Observed")
    z = np.geomspace(x[local].min(), x[local].max(), 100)
    ax.loglog(
        z,
        coefficient * z**expected,
        "k--",
        linewidth=1.3,
        label=f"Slope {expected} reference",
    )
    ax.set_xlabel(r"$\Vert e_n\Vert_2$")
    ax.set_ylabel(r"$\Vert e_{n+1}\Vert_2$")
    ax.set_title(f"TRHiSD-{method} (local fit: {slope:.3f})")
    ax.legend(fontsize=12)
    ax.grid(which="both", alpha=0.2)
    rows.append(
        dict(
            method=f"TRHiSD-{method}",
            iterations=int(r["iterations"]),
            residual=float(np.asarray(r["residual"]).ravel()[-1]),
            error=float(e[-1]),
            local_samples=int(local.sum()),
            loglog_slope=float(slope),
            reference_order=expected,
        )
    )
    with (args.output_dir / f"{method}_error_pairs.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "n",
                "error_n",
                "error_n_plus_1",
                "linear_ratio",
                "quadratic_ratio",
                "local_fit_sample",
            ]
        )
        writer.writerows(zip(range(len(x)), x, y, y / x, y / x**2, local.astype(int)))
with (args.output_dir / "order_summary.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(args.output_dir / f"Rosenbrock_convergence_order.{ext}", dpi=300)
plt.close(fig)
for row in rows:
    print(row)
