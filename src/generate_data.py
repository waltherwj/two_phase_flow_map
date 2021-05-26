"""
This is used to generate the maps of each condition 
and the velocity maps
"""


import numpy as np
from config import Config
from scipy.signal import convolve


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


def detect_edges(parsed_map, kernel=None):
    """get the edge map using the fact that each area is a well defined
    single value section
    """
    if kernel is None:
        # standard kernel
        kernel = np.array(
            [
                [0, 0, -1, 0, 0],
                [0, 0, -1, 0, 0],
                [-1, -1, 8, -1, -1],
                [0, 0, -1, 0, 0],
                [0, 0, -1, 0, 0],
            ]
        )

    edge_map = np.abs(convolve(parsed_map, kernel, mode="same"))

    return edge_map
