""" 
This module dictates how all the maps interact and which ones 
are valid at a certain map location
"""
import numpy as np
from conditions import dispersed_bubbles, stratified
from conditions import bubbly
from conditions import unphysical


def parse_bubbly(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the conditions of the bubbly maps
    """
    # check if bubbly flow can exist
    if bubbly.taylor_velocity_exceeds(
        liquid, gas, pipe
    ) and bubbly.angle_prevents_bubble_migration(liquid, gas, pipe):
        # if it exists calculate it
        gas_void_fraction_bubbly_map = bubbly.gas_void_fraction(
            u_gs, u_ls, liquid, gas, pipe
        )
    else:
        # otherwise the map is just full of False
        gas_void_fraction_bubbly_map = np.full_like(u_gs, False, dtype=bool)

    return gas_void_fraction_bubbly_map


def parse_dispersed_bubble(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the conditions of the dispersed bubble maps
    """
    # get dispersed bubble flow conditions
    gas_void_frac_dispersed_map = dispersed_bubbles.gas_void_fraction(u_gs, u_ls)
    coalescence_map = dispersed_bubbles.bubble_coalescence(
        u_gs, u_ls, liquid, gas, pipe
    )

    # get the correct locations for the gas void fract map
    # based on the coalescence map

    # find the maximum u_gs at which we need to
    # consider both coalescence and void fraction
    maximum_u_gs = u_gs[coalescence_map & gas_void_frac_dispersed_map].max()

    # all locations where dispersed bubbles can exist
    bubble_map = np.full_like(coalescence_map, False)
    bubble_map = gas_void_frac_dispersed_map
    bubble_map[u_gs < maximum_u_gs] = (gas_void_frac_dispersed_map & coalescence_map)[
        u_gs < maximum_u_gs
    ]

    return bubble_map


def parse_stratified(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the conditions of the stratified region
    """
    # get stratified condition region
    stratified_equilibrium_map = stratified.equilibrium_equation(
        u_gs, u_ls, liquid, gas, pipe
    )

    # get the condition of transition into annular
    not_too_steep_map = ~stratified.too_steep_for_stratified(
        u_gs, u_ls, liquid, gas, pipe
    )

    stratified_map = stratified_equilibrium_map & not_too_steep_map
    return stratified_map


def parse_unphysical(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the unphysical condition maps
    """
    impossible_holdup_map = unphysical.total_holdup(u_gs, u_ls, liquid, gas, pipe)

    return impossible_holdup_map


def get_categories_maps(u_gs, u_ls, liquid, gas, pipe):
    """
    calls the other parsing functions to combine all parses into one
    comprehensive map
    """
    bubbly_map = parse_bubbly(u_gs, u_ls, liquid, gas, pipe)
    bubble_map = parse_dispersed_bubble(u_gs, u_ls, liquid, gas, pipe)
    stratified_map = parse_stratified(u_gs, u_ls, liquid, gas, pipe)
    unphysical_map = parse_unphysical(u_gs, u_ls, liquid, gas, pipe)

    # initialize a map with all zeros. Some maps are overlays on the actual map
    category_map = np.full_like(u_ls, np.nan)
    overlay_map = np.full_like(u_ls, np.nan)
    # colors will correspond to these numbers

    # dispersed bubble is true regardless of other conditions
    category_map[bubble_map] = 1

    # stratified if it is stratified
    category_map[stratified_map & np.isnan(category_map)] = 2

    # if bubbly is possible and it is not anything else, it is bubbly (not true but for now)
    category_map[bubbly_map & np.isnan(category_map)] = 3

    # unphysical conditions are always NaN
    # category_map[unphysical_map] = -1
    overlay_map[unphysical_map] = -1

    return category_map, overlay_map
