"""
This module handles the functions that are used to visualize the
maps and other important features
"""
from config import Config
import matplotlib.pyplot as plt
import general
from scipy.ndimage import gaussian_filter
import generate_data


def plot_map(category_map, overlay_map, liquid, gas, pipe, x_ticks=None, y_ticks=None):
    """
    plot a map which contains the representation of a categorical map
    """
    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))

    # get alpha values
    edges = generate_data.detect_edges(category_map)
    alphas = gaussian_filter(edges, sigma=1)
    # scale up to 1
    alphas = alphas / alphas.max()
    # scale min alpha
    alphas[alphas < 0.2] = 0.2

    # get the specific location of the supplied case
    u_gs = general.single_fluid_velocity(gas, pipe)
    u_ls = general.single_fluid_velocity(liquid, pipe)

    c = axs.pcolormesh(
        x_ticks,
        y_ticks,
        category_map,
        shading="gouraud",
        cmap=Config.CMAP,
        alpha=alphas,
    )

    axs.plot(u_gs, u_ls, marker="+", markersize=10)

    fig.colorbar(c)
    axs.set_xticks(x_ticks)
    axs.set_yticks(y_ticks)
    axs.set_xscale("log")
    axs.set_yscale("log")

    return fig, axs


if __name__ == "__main__":
    import numpy as np
    from generate_data import generate_velocity_maps
    from fluids import Liquid, Gas, Pipe
    from parse_maps import get_categories_maps

    total_mass_flow = 0.01
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

    for inclination in [-90]:  # , -80, -30, -1, 0, 1, 30, 85, 90]:
        pipe_temp = Pipe(diameter=5.1e-2, inclination=inclination, roughness=0.00001)
        categories, overlays = get_categories_maps(
            ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
        )

        fig_temp, ax = plot_map(
            categories,
            overlays,
            liq_temp,
            gas_temp,
            pipe_temp,
            x_ticks=ugs_temp[0, :],
            y_ticks=uls_temp[:, 0],
        )

    plt.show(block=True)
