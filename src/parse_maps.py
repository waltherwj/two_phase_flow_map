""" 
This module dictates how all the maps interact and which ones 
are valid at a certain map location
"""
from config import Config
import numpy as np
from conditions import annular, bubbly, dispersed_bubbles, stratified, intermittent


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


def parse_annular(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the conditions of the annular region
    """

    # get liquid stability condition
    stability_map = annular.liquid_stability(u_gs, u_ls, liquid, gas, pipe)
    # get core blockage condition
    core_not_blocked_map = annular.gas_core_blockage(u_gs, u_ls, liquid, gas, pipe)

    annular_map = stability_map & core_not_blocked_map

    return annular_map


def parse_elongated_bubble(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the elongated bubble condition map
    """
    elongated_bubble_map = intermittent.slug_free_of_bubbles(
        u_gs, u_ls, liquid, gas, pipe
    )

    return elongated_bubble_map


def parse_churn(u_gs, u_ls, liquid, gas, pipe):
    """
    parse the slug flow condition condition map
    """
    elongated_bubble_map = intermittent.slug_full_of_bubbles(
        u_gs, u_ls, liquid, gas, pipe
    )

    return elongated_bubble_map


def get_categories_maps(u_gs, u_ls, liquid, gas, pipe):
    """
    calls the other parsing functions to combine all parses into one
    comprehensive map
    """
    bubbly_map = parse_bubbly(u_gs, u_ls, liquid, gas, pipe)
    bubble_map = parse_dispersed_bubble(u_gs, u_ls, liquid, gas, pipe)
    stratified_map = parse_stratified(u_gs, u_ls, liquid, gas, pipe)
    annular_map = parse_annular(u_gs, u_ls, liquid, gas, pipe)
    elongated_bubble_map = parse_elongated_bubble(u_gs, u_ls, liquid, gas, pipe)
    churn_map = parse_churn(u_gs, u_ls, liquid, gas, pipe)

    # initialize a map with all zeros. Some maps are overlays on the actual map
    category_map = np.full_like(u_ls, np.nan)
    overlay_map = np.full_like(u_ls, np.nan)
    # colors will correspond to these numbers

    # dispersed bubble is true regardless of other conditions
    category_map[bubble_map] = Config.CATEGORIES["dispersed bubble"]

    # stratified
    category_map[stratified_map & np.isnan(category_map)] = Config.CATEGORIES[
        "stratified"
    ]

    # annular
    category_map[annular_map & np.isnan(category_map)] = Config.CATEGORIES["annular"]

    # if bubbly is possible and it then it is not an elongated bubble
    if np.any(bubbly_map):
        category_map[bubbly_map & np.isnan(category_map)] = Config.CATEGORIES["bubbly"]
    else:
        # elongated bubble
        category_map[elongated_bubble_map & np.isnan(category_map)] = Config.CATEGORIES[
            "elongated bubble"
        ]

    # slug flow
    category_map[~churn_map & np.isnan(category_map)] = Config.CATEGORIES["slug"]

    # churn flow
    category_map[churn_map & np.isnan(category_map)] = Config.CATEGORIES["churn"]

    return category_map
