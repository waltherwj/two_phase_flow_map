"""
This is used to generate the maps of each condition 
and the velocity maps
"""


import numpy as np
from config import Config


def generate_velocity_maps(
    datapoints=Config.NUMBER_DATAPOINTS,
    min_uls=Config.MIN_ULS,
    max_uls=Config.MAX_ULS,
    min_ugs=Config.MIN_UGS,
    max_ugs=Config.MAX_UGS,
):
    """create maps of corresponding ugs and uls
    where ugs is the x axis and uls is the y axis
    """

    # create the arrays
    ugs_array = np.geomspace(min_ugs, max_ugs, num=datapoints)
    uls_array = np.geomspace(min_uls, max_uls, num=datapoints)

    # tile them up in the correct dimension
    ugs_map = np.tile(ugs_array, (datapoints, 1))
    uls_map = np.tile(uls_array, (datapoints, 1)).T

    return ugs_map, uls_map
