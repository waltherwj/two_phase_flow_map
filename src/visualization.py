"""
This module handles the functions that are used to visualize the
maps and other important features
"""
from scipy.spatial.kdtree import KDTree
from scipy.spatial.qhull import Voronoi
from generate_data import generate_probability_map, refine_velocity_maps
import scipy.interpolate
from scipy.spatial import Delaunay
from conditions import annular
import matplotlib.pyplot as plt
import matplotlib.colors


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
            # norm=matplotlib.colors.LogNorm(thresh),
            shading="nearest",
        )
        fig.colorbar(c)
        axs.set_xticks(x_ticks)
        axs.set_yticks(y_ticks)
        axs.set_xscale("log")
        axs.set_yscale("log")

    return fig, axs


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

    # get new u_ls and u_gs

    # create a mesh with this data
    # x_mesh, y_mesh = np.meshgrid(x_ticks, y_ticks, sparse=True)
    # # interpolate data onto mesh
    # interpolated = scipy.interpolate.griddata(
    #     (x_ticks, y_ticks),
    #     category_map,
    #     (x_mesh, y_mesh),
    #     method="nearest",
    # )

    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))
    data = np.array([x_ticks, y_ticks]).T
    kdtree = KDTree(np.log10(data))
    distances = kdtree.query(np.log10(data), k=3)[0][:, 1:3].sum(axis=1)

    # handle either using axes or not
    c = axs.scatter(
        x_ticks,
        y_ticks,
        c=category_map,
        s=distances * 100,
        alpha=1 - distances / np.max(distances)
        # shading="nearest",
        # cmap="Set1",
        # alpha=0.2,
    )

    axs.set_xscale("log")
    axs.set_yscale("log")

    return fig, axs


if __name__ == "__main__":
    import numpy as np
    from generate_data import generate_velocity_maps, detect_edges, refine_velocity_maps
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

    for inclination in [-90, -80, -30, -1, 0, 1, 30, 85, 90]:
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
        # continuous = detect_edges(categories)
        # continuous, ugs_array, uls_array = generate_probability_map(
        #     continuous, ugs_temp, uls_temp
        # )
        # fig_temp, ax = plot_continuous_symlog(
        #     continuous,
        #     ugs_array,
        #     uls_array,
        # )
        # ax.set_title(f"{inclination}")

        u_gs_refined_temp, u_ls_refined_temp = refine_velocity_maps(
            categories, ugs_temp, uls_temp
        )

        categories, overlays = get_categories_maps(
            u_gs_refined_temp, u_ls_refined_temp, liq_temp, gas_temp, pipe_temp
        )

        plot_categorical_unstructured_map(
            categories.ravel(),
            overlays.ravel(),
            x_ticks=u_gs_refined_temp.ravel(),
            y_ticks=u_ls_refined_temp.ravel(),
        )

        ax.scatter(u_gs_refined_temp, u_ls_refined_temp, s=1)
        # plt.savefig("high_dpis.png")
        plt.show(block=True)
