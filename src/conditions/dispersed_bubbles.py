"""
the functions that define the conditions for flow to be considered
dispersed bubble flow
"""
import numpy as np
import general.friction_factor as friction_factor
import general
import equations.dispersed_bubbles
import fluids


def gas_void_fraction(u_gs, u_ls, critical_value=0.35):
    """condition for gas void fraction
    the realistic gas void fraction is not 0.52,
    as higher gas void fractions than 0.35 are rare
    TODO find citation
    """
    return u_ls > u_gs * (1 - critical_value) / critical_value


def bubble_coalescence(u_gs, u_ls, liquid, gas, pipe):
    """check if the bubbles are large enough to coalesce into
    slugs with consideration for angle
     Barnea 1980
    """

    # local variables for readability
    roughness = pipe.roughness
    sigma = liquid.bubble_surface_tension
    rho_l = liquid.density
    diam = pipe.diameter

    # calculate the lhs
    deformed_bubble_size = equations.dispersed_bubbles.deformed_bubble_critical_size(
        liquid, gas, pipe
    )
    migration_to_top_size = equations.dispersed_bubbles.migration_to_top_critical_size(
        u_gs, u_ls, liquid, gas, pipe
    )
    # make the array
    bubble_crit_diam = np.minimum(migration_to_top_size, deformed_bubble_size)

    # get the mixture values
    mix = fluids.Mix(u_gs, u_ls, liquid, gas, pipe)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)

    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)

    # mixture friction factor
    fric_mix = friction_factor.fang(reynolds_mix, roughness)

    # get the terms for readability
    rhs_1 = 0.725 + 4.15 * np.sqrt(u_gs / u_mix)
    rhs_2 = (sigma / rho_l) ** (3 / 5)
    rhs_3 = ((2 * fric_mix / diam) * (u_mix ** 3)) ** (-2 / 5)
    rhs = rhs_1 * rhs_2 * rhs_3

    # check if the bubbles are few/small enough that they won't coalesce
    return bubble_crit_diam > rhs
