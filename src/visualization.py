"""
This module handles the functions that are used to visualize the
maps and other important features
"""
import matplotlib.pyplot as plt
import matplotlib.colors


def plot_continuous_symlog(array_map, x_ticks=None, y_ticks=None, thresh=0.1):
    """
    plot a map which contains the representation of a
    continuous variable on a symmetrical log scale
    (most scenarios for this type of mapping)
    """
    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))

    # handle either using axes or not
    if (x_ticks is None) or (y_ticks is None):
        axs.pcolormesh(array_map, norm=matplotlib.colors.SymLogNorm(thresh))
    else:
        axs.pcolormesh(
            x_ticks,
            y_ticks,
            array_map,
            norm=matplotlib.colors.SymLogNorm(thresh),
            shading="nearest",
        )
        axs.set_xticks(x_ticks)
        axs.set_yticks(y_ticks)
        axs.set_xscale("log")
        axs.set_yscale("log")

    return fig, axs
