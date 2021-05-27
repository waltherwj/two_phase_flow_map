"""
This is used to generate the maps of each condition 
and the velocity maps
"""


import numpy as np
from config import Config
from scipy import interpolate
from scipy.signal import convolve
from scipy.ndimage import gaussian_filter


def generate_velocity_maps(
    datapoints=Config.NUMBER_REFINED_DATAPOINTS,
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


def generate_velocity_data_close_to_edges(
    rough_map,
    u_gs,
    u_ls,
    refined_datapoints=Config.NUMBER_REFINED_DATAPOINTS,
    sigma=1,
    upscale=10,
):
    """get a new velocity map which takes into account the edges of the
    rough map. The velocity map uses both the original rough map points plus
    a set of refined points which are chosen based on the edges detected on the
    rough.

    This is not used for plotting, but could theoretically be used to have better guesses
    to solve for the transitions
    """
    # detect the edges of the rough map
    edges_map = detect_edges(rough_map)

    # get map with the probability distribution and u_gs and u_ls maps
    # to sample from
    probability_map, u_gs_sample, u_ls_sample = generate_probability_map(
        edges_map, u_gs, u_ls, sigma, upscale
    )
    # tile them up in the correct dimension
    u_gs_map = np.tile(u_gs_sample, (u_gs_sample.size, 1))
    u_ls_map = np.tile(u_ls_sample, (u_ls_sample.size, 1)).T

    # get the original ugs and uls
    u_gs_rough = u_gs.ravel()
    u_ls_rough = u_ls.ravel()

    # get the random sample of indexes using the probability map
    indexes_to_choose = np.arange(probability_map.size)
    indexes = np.random.choice(
        indexes_to_choose,
        size=refined_datapoints ** 2,
        replace=False,
        p=probability_map.ravel(),
    )
    u_gs_refined = u_gs_map.ravel()[indexes]
    u_ls_refined = u_ls_map.ravel()[indexes]

    # join the two
    u_gs_refined = np.append(u_gs_rough, u_gs_refined)
    u_ls_refined = np.append(u_ls_rough, u_ls_refined)

    return u_gs_refined, u_ls_refined


def detect_edges(parsed_map, sigma=0.5):
    """get the edge map using the fact that each area is a well defined
    single value section
    """

    # edges in the x direction
    edges_x = np.diff(parsed_map, axis=0, prepend=parsed_map[:1, :])
    # edges in the y direction
    edges_y = np.diff(parsed_map, axis=1, prepend=parsed_map[:, :1])

    edges_map = ((edges_x != 0) | (edges_y != 0)).astype(float)
    edges_map[edges_map == 0] = np.nan

    return edges_map


def generate_probability_map(edges_map, u_gs, u_ls, sigma=1, upsample=10):
    """generate a probability distribution for points,
    with a high likelihood of putting points close to the
    edges.
    The upsample disctates how fine the probability distribution
    will be. The actual probability distribution is a discrete pdf,
    but the higher the number of points the more it looks like a continuous pdf
    """

    # probability datapoints
    prob_datapoints = Config.NUMBER_REFINED_DATAPOINTS * upsample

    # first smooth out the edges
    probability_map = gaussian_filter(edges_map, sigma=sigma)

    # create an interpolation function with it
    size_x = probability_map.shape[0]
    size_y = probability_map.shape[1]
    x_initial = np.arange(size_x)
    y_initial = np.arange(size_y)
    interpolation_function = interpolate.interp2d(
        x_initial, y_initial, z=probability_map
    )

    # create new vectors with the same range but more points
    x_new = np.linspace(start=0, stop=size_x, num=prob_datapoints)
    y_new = np.linspace(start=0, stop=size_y, num=prob_datapoints)

    # upsample the probability map
    probability_map = interpolation_function(x_new, y_new)

    # upsample the u_gs and u_ls arrays
    # create the arrays
    u_gs_array = np.geomspace(u_gs.min(), u_gs.max(), num=prob_datapoints)
    u_ls_array = np.geomspace(u_ls.min(), u_ls.max(), num=prob_datapoints)

    # normalize it so that it sums to 1
    probability_map /= np.abs(probability_map).sum()

    return probability_map, u_gs_array, u_ls_array
