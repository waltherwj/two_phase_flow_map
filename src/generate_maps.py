"""
This is used to generate the maps of each condition 
and the velocity maps
"""


import numpy as np
from config import Config


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
