"""Shared Plotly colours used by the AutoQY browser interfaces."""


PLOT_BLUE = "#2d6f8e"
PLOT_ORANGE = "#d67b36"
PLOT_PURPLE = "#8a5a9b"
PLOT_BROWN = "#725c42"
PLOT_NEUTRAL = "#6c7280"

ANALYSIS_TRACE_PALETTE = (
    PLOT_BLUE,
    PLOT_ORANGE,
    PLOT_PURPLE,
    PLOT_BROWN,
    PLOT_NEUTRAL,
)

# Okabe–Ito colours; omit yellow, which is faint on the white plot background.
COLORBLIND_TRACE_PALETTE = (
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#D55E00",  # vermilion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#000000",  # black
)
