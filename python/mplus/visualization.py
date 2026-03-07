"""Visualization helpers for Mplus registration results.

All functions return ``matplotlib.figure.Figure`` objects so callers
can ``plt.show()``, ``fig.savefig(…)`` or embed in Jupyter notebooks.

Requires the *visualization* extra::

    pip install mplus-registration[visualization]
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    _HAS_MPL = True
except ImportError:  # pragma: no cover
    _HAS_MPL = False


def _require_mpl() -> None:
    if not _HAS_MPL:
        raise ImportError(
            "matplotlib is required for visualization.  "
            "Install it with: pip install mplus-registration[visualization]"
        )


# ---------------------------------------------------------------------------
# Checkerboard overlay
# ---------------------------------------------------------------------------

def plot_registration_checkerboard(
    fixed: np.ndarray,
    moving_warped: np.ndarray,
    slice_idx: Optional[int] = None,
    axis: int = 2,
    num_squares: int = 8,
    figsize: Tuple[int, int] = (14, 5),
) -> "Figure":
    """Display a checkerboard overlay of fixed and warped-moving images.

    Parameters
    ----------
    fixed : np.ndarray
        3-D fixed image.
    moving_warped : np.ndarray
        3-D warped moving image (same shape as *fixed*).
    slice_idx : int, optional
        Which slice to show along *axis*.  Defaults to the middle slice.
    axis : int
        Axis along which to take the slice (0, 1, or 2).
    num_squares : int
        Number of checker squares per row/column.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _require_mpl()

    if slice_idx is None:
        slice_idx = fixed.shape[axis] // 2

    slc = [slice(None)] * 3
    slc[axis] = slice_idx
    f2d = fixed[tuple(slc)]
    m2d = moving_warped[tuple(slc)]

    # Build checkerboard mask
    rows, cols = f2d.shape
    sq_r = max(rows // num_squares, 1)
    sq_c = max(cols // num_squares, 1)
    rr = (np.arange(rows) // sq_r) % 2
    cc = (np.arange(cols) // sq_c) % 2
    mask = rr[:, None] ^ cc[None, :]

    checker = np.where(mask, f2d, m2d)

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    axes[0].imshow(f2d, cmap="gray")
    axes[0].set_title("Fixed")
    axes[0].axis("off")

    axes[1].imshow(m2d, cmap="gray")
    axes[1].set_title("Warped Moving")
    axes[1].axis("off")

    axes[2].imshow(checker, cmap="gray")
    axes[2].set_title("Checkerboard")
    axes[2].axis("off")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Per-label Dice bar chart
# ---------------------------------------------------------------------------

def plot_label_dice(
    dice_scores: Dict[int, float],
    figsize: Tuple[int, int] = (10, 5),
    color: str = "steelblue",
    title: str = "Per-label Dice Coefficients",
) -> "Figure":
    """Bar chart of per-label Dice coefficients.

    Parameters
    ----------
    dice_scores : dict[int, float]
        Mapping label → Dice ∈ [0, 1].
    figsize, color, title
        Cosmetic parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _require_mpl()

    labels = sorted(dice_scores.keys())
    values = [dice_scores[l] for l in labels]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar([str(l) for l in labels], values, color=color)
    ax.set_xlabel("Label")
    ax.set_ylabel("Dice Coefficient")
    ax.set_title(title)
    ax.set_ylim(0, 1.05)
    ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.5)

    # Annotate each bar with the value
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{val:.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Metric convergence plot
# ---------------------------------------------------------------------------

def plot_metric_convergence(
    metric_history: Union[List[float], Dict[str, List[float]]],
    figsize: Tuple[int, int] = (10, 5),
    title: str = "Metric Convergence",
) -> "Figure":
    """Line plot of metric value(s) over optimizer iterations.

    Parameters
    ----------
    metric_history : list[float] or dict[str, list[float]]
        If a list, the combined metric per iteration.
        If a dict, a history per sub-metric name (MI, NGF, …).
    figsize, title
        Cosmetic parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _require_mpl()

    fig, ax = plt.subplots(figsize=figsize)

    if isinstance(metric_history, dict):
        for name, vals in metric_history.items():
            ax.plot(vals, label=name)
        ax.legend()
    else:
        ax.plot(metric_history, label="Combined", color="steelblue")

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Metric Value")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Slice comparison (before / after)
# ---------------------------------------------------------------------------

def plot_before_after(
    fixed: np.ndarray,
    moving: np.ndarray,
    warped: np.ndarray,
    slice_idx: Optional[int] = None,
    axis: int = 2,
    figsize: Tuple[int, int] = (18, 5),
) -> "Figure":
    """Side-by-side comparison: fixed, original moving, warped.

    Parameters
    ----------
    fixed, moving, warped : np.ndarray
        3-D images.
    slice_idx : int, optional
        Slice along *axis*.  Defaults to middle.
    axis : int
        Slicing axis.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _require_mpl()

    if slice_idx is None:
        slice_idx = fixed.shape[axis] // 2

    slc = [slice(None)] * 3
    slc[axis] = slice_idx

    images = [fixed[tuple(slc)], moving[tuple(slc)], warped[tuple(slc)]]
    titles = ["Fixed", "Moving (original)", "Moving (warped)"]

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    for ax, img, t in zip(axes, images, titles):
        ax.imshow(img, cmap="gray")
        ax.set_title(t)
        ax.axis("off")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Difference map
# ---------------------------------------------------------------------------

def plot_difference_map(
    fixed: np.ndarray,
    warped: np.ndarray,
    slice_idx: Optional[int] = None,
    axis: int = 2,
    figsize: Tuple[int, int] = (8, 6),
    cmap: str = "RdBu_r",
) -> "Figure":
    """Heat-map of the intensity difference (fixed − warped).

    Parameters
    ----------
    fixed, warped : np.ndarray
        3-D images.
    slice_idx, axis, figsize, cmap
        Display parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _require_mpl()

    if slice_idx is None:
        slice_idx = fixed.shape[axis] // 2

    slc = [slice(None)] * 3
    slc[axis] = slice_idx
    diff = fixed[tuple(slc)].astype(np.float64) - warped[tuple(slc)].astype(np.float64)

    vmax = max(abs(diff.min()), abs(diff.max())) or 1.0

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(diff, cmap=cmap, vmin=-vmax, vmax=vmax)
    fig.colorbar(im, ax=ax, label="Fixed − Warped")
    ax.set_title("Difference Map")
    ax.axis("off")
    fig.tight_layout()
    return fig
