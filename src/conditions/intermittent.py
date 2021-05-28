"""
the functions that define the conditions for flow to be considered
intermittent flow
"""

import equations


def slug_free_of_bubbles(u_gs, u_ls, liquid, gas, pipe):
    """condition for if liquid slug is free of entrained bubbles"""

    # calculate the holdup of gas inside of the liquid slug
    gas_holdup_in_slug = equations.intermittent.liquid_slug_gas_holdup(
        u_gs, u_ls, liquid, gas, pipe
    )
    liquid_holdup = 1 - gas_holdup_in_slug

    # the condition
    return liquid_holdup >= 1


def slug_full_of_bubbles(u_gs, u_ls, liquid, gas, pipe):
    """condition for if liquid slug is at the maximum packing of
    of entrained bubbles where the slug collapses
    """
    # calculate the holdup of gas inside of the liquid slug
    gas_holdup_in_slug = equations.intermittent.liquid_slug_gas_holdup(
        u_gs, u_ls, liquid, gas, pipe
    )
    liquid_holdup = 1 - gas_holdup_in_slug

    # the condition
    return liquid_holdup <= 0.48
