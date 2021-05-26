"""
This is used to generate the maps of each condition 
and the velocity maps
"""


import numpy as np
from numpy.core.numeric import full_like
from config import Config
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter


def generate_velocity_maps(
    datapoints=Config.NUMBER_DATAPOINTS,
    min_u_ls=Config.MIN_ULS,
    max_u_ls=Config.MAX_ULS,
    min_u_gs=Config.MIN_UGS,
    max_u_gs=Config.MAX_UGS,
):
    """create maps of corresponding u_gs and u_ls
    where u_gs is the x axis and u_ls is the y axis
    """

    # create the arrays
    u_gs_array = np.geomspace(min_u_gs, max_u_gs, num=datapoints)
    u_ls_array = np.geomspace(min_u_ls, max_u_ls, num=datapoints)

    # tile them up in the correct dimension
    u_gs_map = np.tile(u_gs_array, (datapoints, 1))
    u_ls_map = np.tile(u_ls_array, (datapoints, 1)).T

    return u_gs_map, u_ls_map


def refine_velocity_maps(rough_map, u_gs, u_ls):
    """get a new velocity map which takes into account the edges of the
    rough map
    """
    pass


def detect_edges(parsed_map):
    """get the edge map using the fact that each area is a well defined
    single value section
    """

    # edges in the x direction
    edges_x = np.diff(parsed_map, axis=0, prepend=parsed_map[:1, :])
    # edges in the y direction
    edges_y = np.diff(parsed_map, axis=1, prepend=parsed_map[:, :1])

    edges_map = ((edges_x != 0) | (edges_y != 0)).astype(float)
    # edge_map = np.abs(convolve(parsed_map, kernel, mode="same"))
    # edge_map[edge_map < 1] = 1

    return edges_map


def generate_probability_map(edges_map):
    """generate a probability distribution for points
    that has a high likelihood of putting points close to the
    edges
    """

    # first smooth out the edges
    probability_map = gaussian_filter(edges_map, sigma=3)

    return probability_map
