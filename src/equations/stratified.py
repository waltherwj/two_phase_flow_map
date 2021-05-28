"""
the functions that are used to calculate the conditions for flow to be considered
intermittent flow
"""
import numpy as np
from general.general import Geometry


def wave_growth(crit_height, u_gs, liquid, gas, pipe):
    """f(x) = 0 formulation to find the critical fluid height
    at which waves will grow
    Taitel&Duckler 1976
    """

    # get the modified froude number
    froude = modified_froude(u_gs, liquid, gas, pipe)

    # fix broken values
    crit_height[crit_height > 1] = 1
    crit_height[crit_height < 0] = 0

    # non dimensional values
    tilde = Geometry(crit_height, non_dimensional=True)

    # derivative
    dadh_l_tilde = ((2 - 2 * crit_height) * crit_height) / np.sqrt(
        (1 - crit_height) * crit_height
    )

    # divide lhs in terms for readability
    lhs_1 = 1 / ((1 - crit_height) ** 2)
    lhs_2 = (tilde.vel_g ** 2) * dadh_l_tilde / tilde.area_g
    lhs = (froude ** 2) * lhs_1 * lhs_2

    return lhs - 1


def modified_froude(u_gs, liquid, gas, pipe):
    """calculate the modified froude modified by density ratio"""

    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density
    diam = pipe.diameter
    grav = pipe.gravity
    beta = pipe.inclination

    # terms for readability
    froude_1 = np.sqrt(rho_g / (rho_l - rho_g))
    froude_2 = u_gs / np.sqrt(diam * grav * np.cos(beta))
    froude = froude_1 * froude_2

    return froude
