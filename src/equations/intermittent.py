"""
the functions that are used to calculate the conditions for flow to be considered
intermittent flow
"""
import numpy as np
from . import dispersed_bubbles
from general import friction_factor
import general
import fluids


def liquid_slug_gas_holdup(u_gs, u_ls, liquid, gas, pipe):
    """
    calculate the slug holdup based on mixed properties
    Barnea 1987
    """
    # local variables for readability
    rho_l = liquid.density
    sigma = liquid.bubble_surface_tension
    diam = pipe.diameter

    # get the mixture values
    mix = fluids.Mix(u_gs, u_ls, liquid, gas, pipe)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)
    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)
    # mixture friction factor
    fric_mix = friction_factor.niazkar_and_churchill(reynolds_mix, pipe.roughness)

    # get the critical size
    critical_diam = dispersed_bubbles.deformed_bubble_critical_size(liquid, gas, pipe)

    # the expression split in terms for readability
    term_1 = critical_diam * (2 * fric_mix * (u_mix ** 3) / diam) ** (2 / 5)
    term_2 = (rho_l / sigma) ** (3 / 5)

    holdup = 0.058 * (term_1 * term_2 - 0.725) ** 2

    # recuperate the sign that was lost in the maths
    sign = np.sign(term_1 * term_2 - 0.725)

    # negative holdups do not make sense, but this equation is used in ways
    # that requires the positive and negative numbers to be available
    holdup = sign * holdup

    return holdup
