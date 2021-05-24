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

    # all locations where dispersed bubbles can exist
    bubble_map = coalescence_map & (gas_void_frac_dispersed_map)

    return bubble_map


def parse_stratified(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the conditions of the stratified region
    """
    # get stratified condition region
    stratified_equilibrium_map = stratified.equilibrium_equation(
        u_gs, u_ls, liquid, gas, pipe
    )

    stratified_map = stratified_equilibrium_map
    return stratified_map


def parse_unphysical(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the unphysical condition maps
    """
    impossible_holdup_map = unphysical.total_holdup(u_gs, u_ls, liquid, gas, pipe)

    return impossible_holdup_map


def get_categories_map(u_gs, u_ls, liquid, gas, pipe):
    """
    calls the other parsing functions to combine all parses into one
    comprehensive map
    """
    bubbly_map = parse_bubbly(u_gs, u_ls, liquid, gas, pipe)
    bubble_map = parse_dispersed_bubble(u_gs, u_ls, liquid, gas, pipe)
    stratified_map = parse_stratified(u_gs, u_ls, liquid, gas, pipe)
    unphysical_map = parse_unphysical(u_gs, u_ls, liquid, gas, pipe)

    # initialize a map with all zeros
    category_map = np.zeros_like(unphysical_map, dtype=np.float32)
    # colors will correspond to these numbers

    # if bubbly is true, dispersed bubble isn't
    if np.any(bubbly_map):
        bubble_map = bubbly_map
        category_map[bubbly_map] = 2
    else:
        category_map[bubble_map] = 1

    # stratified superseeds bubbles according to the maps
    # barnea's paper algorithm sees to say the opposite,
    # but maybe they were just bad at algorithms back then
    category_map[stratified_map] = 3

    # unphysical conditions are always NaN
    # category_map[unphysical_map] = np.nan

    return category_map
