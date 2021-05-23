"""
This module handles the functions that are used to visualize the
maps and other important features
"""
import matplotlib.pyplot as plt
import matplotlib.colors

from conditions import dispersed_bubbles
import general


def plot_continuous_symlog(array_map, x_ticks=None, y_ticks=None, thresh=1e-5):
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
        c = axs.pcolormesh(
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

    return fig, c


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
    import conditions.dispersed_bubbles

    total_mass_flow = 0.001
    quality = 0.01
    liq_massflow = total_mass_flow * (1 - quality)
    gas_massflow = total_mass_flow * quality

    ugs_temp, uls_temp = generate_velocity_maps()
    liq_temp = Liquid(
        density=998,
        bubble_surface_tension=0.073,
        mass_flowrate=liq_massflow,
        dynamic_viscosity=8.9e-4,
    )
    gas_temp = Gas(density=1.225, mass_flowrate=gas_massflow, dynamic_viscosity=18.3e-6)

    pipe_temp = Pipe(diameter=5.1e-2, inclination=0, roughness=0.001)
    void_frac_map = general.fluid_area_ratio(ugs_temp, gas_temp, pipe_temp)

    fig_temp, ax = plot_continuous_symlog(
        void_frac_map,
        x_ticks=ugs_temp[0, :],
        y_ticks=uls_temp[:, 0],
    )
    fig_temp.colorbar(ax)

    # for inclination in [-90, -80, -30, -1, 0, 1, 30, 80, 90]:
    #     pipe_temp = Pipe(diameter=5.1e-2, inclination=inclination, roughness=0.001)
    #     gas_void_fraction_map = conditions.dispersed_bubbles.gas_void_fraction(
    #         ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
    #     )
    #     coalescence_map = conditions.dispersed_bubbles.bubble_coalescence(
    #         ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
    #     )
    #     fig_temp, ax = plot_categorical_map(
    #         coalescence_map & gas_void_fraction_map,
    #         x_ticks=ugs_temp[0, :],
    #         y_ticks=uls_temp[:, 0],
    #     )
    #     ax.set_title(inclination)
    plt.show(block=True)
