"""
the functions that define the conditions for flow to be considered
intermittent flow
"""
import numpy as np
from scipy.optimize import newton
from fluids import Mix
import general
from general import friction_factor
from .dispersed_bubbles import deformed_bubble_critical_size


def slug_free_of_bubbles(u_gs, u_ls, liquid, gas, pipe):
    """condition for if liquid slug is free of entrained bubbles"""

    # calculate the holdup of gas inside of the liquid slug
    # get the mixture values
    gas_holdup_in_slug = liquid_slug_gas_holdup(u_gs, u_ls, liquid, gas, pipe)
    liquid_holdup = 1 - gas_holdup_in_slug

    # the condition
    return liquid_holdup >= 1


def slug_full_of_bubbles(u_gs, u_ls, liquid, gas, pipe):
    """condition for if liquid slug is at the maximum packing of
    of entrained bubbles where the slug collapses
    """
    # calculate the holdup of gas inside of the liquid slug
    initial_holdup = liquid_slug_gas_holdup(u_gs, u_ls, liquid, gas, pipe)
    gas_holdup_in_slug = newton(
        solve_for_slug_holdup, initial_holdup, args=(u_gs, u_ls, liquid, gas, pipe)
    )
    liquid_holdup = 1 - gas_holdup_in_slug

    liquid_holdup = 1 - gas_holdup_in_slug

    # the condition
    return liquid_holdup <= 0.48


def solve_for_slug_holdup(alpha_g, u_gs, u_ls, liquid, gas, pipe):
    """f(x)= 0  formulation to find the critical slug holdup everywhere"""

    # local variables for readability
    rho_l = liquid.density
    sigma = liquid.bubble_surface_tension
    diam = pipe.diameter

    # get the mixture values
    mix = Mix(u_gs, u_ls, liquid, gas, pipe, gas_holdup=alpha_g)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)
    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)
    # mixture friction factor
    fric_mix = friction_factor.niazkar_and_churchill(reynolds_mix, pipe.roughness)
    # get the critical size
    critical_diam = deformed_bubble_critical_size(liquid, gas, pipe)

    # the expression split in terms for readability
    term_1 = critical_diam * (2 * fric_mix * (u_mix ** 3) / diam) ** (2 / 5)
    term_2 = (rho_l / sigma) ** (3 / 5)

    holdup = 0.058 * (term_1 * term_2 - 0.725) ** 2

    return alpha_g - holdup


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
    mix = Mix(u_gs, u_ls, liquid, gas, pipe)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)
    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)
    # mixture friction factor
    fric_mix = friction_factor.niazkar_and_churchill(reynolds_mix, pipe.roughness)

    # get the critical size
    critical_diam = deformed_bubble_critical_size(liquid, gas, pipe)

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
