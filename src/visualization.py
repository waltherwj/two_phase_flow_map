"""
This module handles the functions that are used to visualize the
maps and other important features
"""
from config import Config
from scipy.spatial.kdtree import KDTree
from scipy.spatial.qhull import Voronoi
from generate_data import generate_probability_map
from scipy.spatial import Delaunay
from conditions import annular
import matplotlib.pyplot as plt
import matplotlib.colors


def plot_continuous_symlog(
    array_map, x_ticks=None, y_ticks=None, axs=None, thresh=1e-5
):
    """
    plot a map which contains the representation of a
    continuous variable on a symmetrical log scale
    (most scenarios for this type of mapping)
    """
    # initialize the figure
    if axs is None:
        fig, axs = plt.subplots(figsize=(7, 7))

    # handle either using axes or not
    if (x_ticks is None) or (y_ticks is None):
        axs.pcolormesh(array_map, norm=matplotlib.colors.SymLogNorm(thresh))
    else:
        c = axs.pcolormesh(
            x_ticks,
            y_ticks,
            array_map,
            shading="gouraud",
            cmap="binary",
            antialiased=True,
        )
        if axs is None:
            fig.colorbar(c)
        axs.set_xticks(x_ticks)
        axs.set_yticks(y_ticks)
        axs.set_xscale("log")
        axs.set_yscale("log")

    return axs


def plot_categorical_map(category_map, overlay_map, x_ticks=None, y_ticks=None):
    """
    plot a map which contains the representation of a categorical map
    """
    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))

    # handle either using axes or not
    if (x_ticks is None) or (y_ticks is None):
        axs.pcolormesh(category_map)
    else:
        c = axs.pcolormesh(
            x_ticks, y_ticks, category_map, shading="gouraud", cmap="Set1", alpha=0.9
        )
        fig.colorbar(c)
        axs.pcolormesh(
            x_ticks, y_ticks, overlay_map, shading="nearest", cmap="Dark2", alpha=0.5
        )
        axs.set_xticks(x_ticks)
        axs.set_yticks(y_ticks)
        axs.set_xscale("log")
        axs.set_yscale("log")

    return fig, axs


def plot_categorical_unstructured_map(category_map, overlay_map, x_ticks, y_ticks):
    """
    plot a map which contains the representation of a categorical map
    when the data is unstructured
    """
    # TODO implement mapping nearest color to background mesh
    # get new u_ls and u_gs
    u_ls_res = np.geomspace(y_ticks.min(), y_ticks.max(), num=Config.RESOLUTION)
    u_gs_res = np.geomspace(x_ticks.min(), x_ticks.max(), num=Config.RESOLUTION)
    # tile them up in the correct dimension and get an array with all location combinations
    u_gs_map = np.tile(u_gs_res, (Config.RESOLUTION, 1)).ravel()
    u_ls_map = np.tile(u_ls_res, (Config.RESOLUTION, 1)).T.ravel()
    data_map = np.array([u_gs_map, u_ls_map]).T

    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))
    data = np.array([x_ticks, y_ticks]).T
    kdtree = KDTree(np.log10(data))
    distances, indices = kdtree.query(np.log10(data_map), k=2)
    # get rid of the zero distance (identity) indexes
    indices = indices[:, 1:]
    # get the array of colors
    colors_array = np.full(Config.RESOLUTION * Config.RESOLUTION, np.nan)
    # fill the array
    colors_array = category_map[indices]
    # shape it into the correct dimension
    colors_array = colors_array.reshape((Config.RESOLUTION, Config.RESOLUTION))
    colors_array = np.rot90(colors_array.T, k=4)
    # plot the image
    axs.hexbin(u_gs_map, u_ls_res, colors_array)

    # handle either using axes or not
    # c = axs.scatter(
    #     x_ticks,
    #     y_ticks,
    #     c=category_map
    #     # shading="nearest",
    #     # cmap="Set1",
    #     # alpha=0.2,
    # )

    # axs.set_xscale("log")
    # axs.set_yscale("log")

    return fig, axs


if __name__ == "__main__":
    import numpy as np
    from generate_data import (
        generate_velocity_maps,
        detect_edges,
        generate_velocity_data_close_to_edges,
    )
    from fluids import Liquid, Gas, Pipe
    from parse_maps import get_categories_maps
    from conditions import stratified, dispersed_bubbles

    total_mass_flow = 0.001
    quality = 0.001
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

    for inclination in [80, -80]:  # [-90, -80, -30, -1, 0, 1, 30, 85, 90]:
        pipe_temp = Pipe(diameter=5.1e-2, inclination=inclination, roughness=0.00001)
        categories, overlays = get_categories_maps(
            ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
        )

        fig_temp, ax = plot_categorical_map(
            categories,
            overlays,
            x_ticks=ugs_temp[0, :],
            y_ticks=uls_temp[:, 0],
        )
        continuous = detect_edges(categories)
        ax = plot_continuous_symlog(continuous, ugs_temp[0, :], uls_temp[:, 0], axs=ax)
        ax.set_title(f"{inclination}")

        # ax.scatter(u_gs_refined_temp, u_ls_refined_temp, s=1)
        # plt.savefig("high_dpis.png")
    plt.show(block=True)
