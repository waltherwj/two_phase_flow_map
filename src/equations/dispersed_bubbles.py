"""
the functions that are used to calculate the conditions for flow to be considered
bubbly flow

"""
import numpy as np
import general
from general import friction_factor
import fluids


def deformed_bubble_critical_size(liquid, gas, pipe):
    """calculate the critical size of a bubble for it to become deformed
    eq. 6 Barnea 1987
    """
    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    grav = pipe.gravity

    # critical bubble size
    diam_crit_deformed = 2 * np.sqrt(0.4 * sigma / ((rho_l - rho_g) * grav))

    return diam_crit_deformed


def migration_to_top_critical_size(u_gs, u_ls, liquid, gas, pipe):
    """calculate the critical size below which a bubble can't travel
    to the upper part of pipe
    eq. 6 Barnea 1987
    """
    # get the mixture
    mix = fluids.Mix(u_gs, u_ls, liquid, gas, pipe)

    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density
    roughness = pipe.roughness
    grav = pipe.gravity
    beta = pipe.inclination

    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)

    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)

    # mixture friction factor
    friction_mix = friction_factor.fang(reynolds_mix, roughness)

    # critical bubble size
    diam_crit_migration = (
        (3 / 8) * (rho_l / (rho_l - rho_g)) * (friction_mix * (u_mix ** 2))
    ) / (grav * np.cos(beta))

    return diam_crit_migration
