"""AGU-style figure settings shared by the plotting scripts.

Follows the AGU publication figure guidelines: Arial/Helvetica sans-serif
text, >= 8 pt font at final size, >= 300 dpi raster output, TrueType fonts
embedded in vector output, and lowercase panel labels "(a)", "(b)", ...
"""

from __future__ import annotations

import matplotlib

AGU_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 9,
    "figure.titlesize": 12,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # embed TrueType so text stays editable
    "ps.fonttype": 42,
}


def apply_agu_style() -> None:
    matplotlib.rcParams.update(AGU_RC)


def panel_label(ax, letter: str, x: float = 0.0, y: float = 1.0,
                outside: bool = False) -> None:
    """Bold '(a)'-style label at the top-left corner of the axes box.

    ``outside=True`` sets it just above the frame instead of inside it.
    """
    if outside and y == 1.0:
        y = 1.01
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="left", va="bottom" if outside else "top",
            clip_on=False)
