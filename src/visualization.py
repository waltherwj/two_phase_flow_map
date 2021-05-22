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


def plot_categorical_map(array_map, x_ticks=None, y_ticks=None):
    """
    plot a map which contains the representation of a categorical map
    """
    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))

    # handle either using axes or not
    if (x_ticks is None) or (y_ticks is None):
        axs.pcolormesh(array_map)
    else:
        axs.pcolormesh(
            x_ticks,
            y_ticks,
            array_map,
            shading="nearest",
        )
        axs.set_xticks(x_ticks)
        axs.set_yticks(y_ticks)
        axs.set_xscale("log")
        axs.set_yscale("log")

    return fig, axs


if __name__ == "__main__":
    import numpy as np
    from generate_maps import generate_velocity_maps
    from fluids import Liquid, Gas, Pipe
    import conditions.slug_flow

    ugs_temp, uls_temp = generate_velocity_maps()
    liq_temp = Liquid(density=998, bubble_surface_tension=0.073)
    gas_temp = Gas(density=1.225)
    for inclination in np.linspace(-90, 90, num=5):
        pipe_temp = Pipe(diameter=5.1e-2, inclination=inclination)
        gas_void_fraction_map = conditions.slug_flow.gas_void_fraction(
            ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
        )
        fig_temp, ax = plot_categorical_map(
            gas_void_fraction_map, x_ticks=ugs_temp[0, :], y_ticks=uls_temp[:, 0]
        )
        ax.set_title(inclination)
    plt.show(block=True)
