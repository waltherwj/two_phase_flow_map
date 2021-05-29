"""
This module handles the functions that are used to visualize the
maps and other important features
"""
from config import Config
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import general
from scipy.ndimage import gaussian_filter
import generate_data
import numpy as np


def plot_map(category_map, liquid, gas, pipe, u_gs_map, u_ls_map):
    """
    plot a map which contains the representation of a categorical map
    """
    # initialize the figure
    fig, axs = plt.subplots(figsize=(7, 7))

    # get the axis values
    x_ticks = u_gs_map[0, :]
    y_ticks = u_ls_map[:, 0]

    # get alpha values
    edges = generate_data.detect_edges(category_map)
    alphas = gaussian_filter(edges, sigma=1)
    # scale up to 1
    alphas = alphas / alphas.max()
    # scale min alpha
    min_alpha = 0.2
    alphas[alphas < min_alpha] = min_alpha

    # get the specific location of the supplied case
    u_gs = general.single_fluid_velocity(gas, pipe)
    u_ls = general.single_fluid_velocity(liquid, pipe)

    # get the legend
    legend_elements = []
    for category, value in Config.CATEGORIES.items():
        if np.any(category_map == value):
            color = list(Config.CMAP(value))
            color[-1] = min(min_alpha * 3, 1)
            legend_elements.append(Patch(facecolor=color, label=category))

    axs.pcolormesh(
        x_ticks,
        y_ticks,
        category_map,
        shading="gouraud",
        cmap=Config.CMAP,
        alpha=alphas,
    )

    # the current value marker
    axs.plot(u_gs, u_ls, marker="x", color="black", mew=2, markersize=10)
    legend_elements.append(
        mlines.Line2D(
            [],
            [],
            color="black",
            marker="x",
            mew=2,
            linestyle="None",
            markersize=10,
            label="This scenario",
        )
    )

    # plot the legend
    axs.legend(handles=legend_elements, bbox_to_anchor=(1.04, 1), loc="upper left")

    # set the ticks, label, title and scale
    axs.set_xticks(x_ticks)
    axs.set_yticks(y_ticks)
    axs.set_xscale("log")
    axs.set_yscale("log")
    axs.set_xlabel(r"$U_{Gs}$")
    axs.set_ylabel(r"$U_{Ls}$")
    axs.set_title(f"Inclination = ${pipe.inclination*180/np.pi:.1f}{{\degree}}$")
    plt.tight_layout()

    return fig, axs
