"""Plot four methods, retaining only the BB run capped at 10*beta0."""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (
    LogLocator,
    LogFormatterMathtext,
    NullFormatter,
    ScalarFormatter,
)
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(
    description="Plot archived or recomputed Rosenbrock results."
)
parser.add_argument("--data-dir", type=Path, default=HERE / "data")
parser.add_argument("--output-dir", type=Path, default=HERE / "figures")
args = parser.parse_args()
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
saved = loadmat(args.data_dir / "d2_all_methods.mat", simplify_cells=True)
archived = loadmat(args.data_dir / "plot_data.mat", simplify_cells=True)
assert saved["runs"][1]["options"]["bbmax"] == 0.01
left = archived["D2"]["lines"]
right = archived["D1000"]["lines"]
BB_COLOR = left[1]["Color"]

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.size": 10,
        "mathtext.fontset": "stix",
        "svg.fonttype": "path",
        "axes.linewidth": 0.65,
        "axes.labelsize": 20,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)
fig = plt.figure(figsize=(11.92, 5.77), dpi=100)
positions = [
    [0.080, 0.155, 0.375, 0.660],
    [0.600, 0.155, 0.375, 0.660],
    [0.275, 0.685, 0.160, 0.110],
    [0.795, 0.685, 0.160, 0.110],
]
axes = [fig.add_axes(pos) for pos in positions]


def panel(ax, data, is_left, inset=False):
    handles = []
    for j, curve in enumerate(data):
        if is_left:
            run = saved["runs"][j]
            x = np.arange(int(run["iterations"]) + 1)
            y = np.asarray(run["error"]).ravel()
        else:
            x = np.asarray(curve["X"]).ravel()
            y = np.asarray(curve["Y"]).ravel()
        name = "HiSD-BB" if j == 1 else curve["Name"]
        (line,) = ax.plot(
            x,
            y,
            color=BB_COLOR if j == 1 else curve["Color"],
            linewidth=1,
            linestyle="-" if j == 1 else curve["LineStyle"],
            label=name,
        )
        assert np.array_equal(line.get_xdata(), x)
        assert np.array_equal(line.get_ydata(), y)
        handles.append(line)
    ax.set_yscale("log")
    ax.set_ylim([1e-6, 1] if is_left else archived["D1000"]["YLim"])
    if inset and not is_left:
        # Fit the positive errors over the requested first 1000 updates.
        visible = np.concatenate(
            [
                line.get_ydata()[(line.get_xdata() >= 0) & (line.get_xdata() <= 1000)]
                for line in handles
            ]
        )
        ax.set_ylim(visible.min() / 1.4, visible.max() * 1.3)
    ax.set_xlim(
        ([0, 40] if is_left else [-50, 1000])
        if inset
        else ([-5, 310] if is_left else archived["D1000"]["XLim"])
    )
    if inset:
        ax.set_xticks(np.arange(0, 41, 10) if is_left else np.arange(0, 1001, 250))
        ax.set_yticks([1e-6, 1e-3, 1] if is_left else [1e-2, 1e-1, 1])
    else:
        if is_left:
            ax.set_xticks(np.arange(0, 301, 100))
        if not is_left:
            ax.set_xticks(np.arange(0, 100001, 20000))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((4, 4))
            ax.xaxis.set_major_formatter(formatter)
        ax.set_yticks(10.0 ** np.arange(-6 if is_left else -4, 1, 2))
    ax.yaxis.set_major_formatter(LogFormatterMathtext())
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(True, which="major", color=".78", alpha=0.45, linewidth=0.5)
    ax.grid(True, which="minor", color=".82", alpha=0.30, linewidth=0.4, linestyle=":")
    ax.tick_params(
        which="both",
        direction="in",
        top=True,
        right=True,
        labelsize=11 if inset else 17,
        width=0.55,
        length=2 if inset else 3,
    )
    for spine in ax.spines.values():
        spine.set_color(".35")
        spine.set_linewidth(0.65)
    label_size = 12 if inset else 20
    ax.set_xlabel(
        r"Iteration number $n$", fontsize=label_size, labelpad=2 if inset else 4
    )
    ax.set_ylabel(
        r"$\Vert x^n-x^*\Vert_2$", fontsize=label_size, labelpad=2 if inset else 4
    )
    ax.yaxis.set_label_coords(-0.26 if inset else -0.120, 0.5)
    return handles


handles = panel(axes[0], left, True)
panel(axes[1], right, False)
panel(axes[2], left, True, True)
panel(axes[3], right, False, True)
for ax, title in zip(axes[:2], [r"(a) $d=2$", r"(b) $d=1000$"]):
    ax.set_title(title, fontsize=16, fontweight="normal", pad=10)

order = [0, 1, 2, 3]
legend = fig.legend(
    [handles[i] for i in order],
    [handles[i].get_label() for i in order],
    loc="upper center",
    bbox_to_anchor=(0.51, 0.970),
    ncol=4,
    prop={"family": "Arial", "size": 14},
    frameon=True,
    fancybox=False,
    edgecolor=".4",
    facecolor="white",
    borderpad=0.25,
    handlelength=1.8,
    handletextpad=0.35,
    columnspacing=0.9,
)

# The same style object is applied to both main axes and to both insets.
for a, b in [(axes[0], axes[1]), (axes[2], axes[3])]:
    assert a.xaxis.label.get_fontsize() == b.xaxis.label.get_fontsize()
    assert a.yaxis.label.get_fontsize() == b.yaxis.label.get_fontsize()
    assert a.xaxis.label.get_fontfamily() == b.xaxis.label.get_fontfamily()
    assert a.spines["left"].get_linewidth() == b.spines["left"].get_linewidth()

assert axes[0].get_ylim() == axes[2].get_ylim() == (1e-6, 1.0)
assert axes[2].get_xlim() == (0.0, 40.0)
assert axes[3].get_xlim() == (-50.0, 1000.0)
assert all(len(ax.lines) == 4 for ax in axes)
fig.canvas.draw()
# Verify the inset including its labels does not hide any main-curve sample.
renderer = fig.canvas.get_renderer()
legend_bounds = legend.get_window_extent(renderer)
assert legend_bounds.x0 >= 0 and legend_bounds.x1 <= fig.bbox.width
for ax in axes[:2]:
    assert legend_bounds.y0 > ax.title.get_window_extent(renderer).y1
for ax in axes:
    for label in [ax.xaxis.label, ax.yaxis.label]:
        bounds = label.get_window_extent(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
for base, inset in [(axes[0], axes[2]), (axes[1], axes[3])]:
    box = inset.get_tightbbox(renderer)
    # Include tick text and axis labels when checking containment.
    main_box = base.get_window_extent(renderer)
    assert box.x0 >= main_box.x0 and box.x1 <= main_box.x1
    assert box.y0 >= main_box.y0 and box.y1 <= main_box.y1
    for line in base.lines:
        xy = base.transData.transform(
            np.column_stack([line.get_xdata(), line.get_ydata()])
        )
        inside = (
            (xy[:, 0] >= box.x0)
            & (xy[:, 0] <= box.x1)
            & (xy[:, 1] >= box.y0)
            & (xy[:, 1] <= box.y1)
        )
        assert (
            not inside.any()
        ), f"Inset labels cover {line.get_label()}: bbox={box.bounds}, x={line.get_xdata()[inside][:3]}"

fig.savefig(OUT / "Rosenbrock_BB_cap_comparison.png", dpi=300)
fig.savefig(OUT / "Rosenbrock_BB_cap_comparison.svg")
fig.savefig(OUT / "Rosenbrock_BB_cap_comparison.pdf")
plt.close(fig)
print("Only cap10 BB retained in all four panels; legend label: HiSD-BB.")
print(
    "Right inset: first 1000 updates, log-y limits fitted to visible errors; main ranges preserved."
)
