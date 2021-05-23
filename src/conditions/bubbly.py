"""
the functions that define the conditions for flow to be considered
bubbly flow
"""

import numpy as np


def taylor_velocity_exceeds(liquid, gas, pipe):
    """check if taylor bubble velocity exceeds bubble velocity
    which is necessary for bubbly flow to exist
        Taitel et al. 1980
    """
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    grav = pipe.gravity
    diam = pipe.diameter

    # calculate the right hand side
    rhs = np.sqrt(19 * (rho_l - rho_g) * sigma / ((rho_l ** 2) * grav))

    # return the conditon
    return diam > rhs


def angle_prevents_bubble_migration(liquid, gas, pipe, lift_coef=None, gamma=None):
    """check if the angle of inclination prevents bubble migration
    to the top wall of the pipe
        Barnea et al. 1985
    """
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    grav = pipe.gravity
    beta = pipe.inclination
    diam = pipe.diameter

    # constants
    # lif coefficient of bubble
    if lift_coef is None:
        # average value from Streeter 1961
        lift_coef = 0.8
    if gamma is None:
        # 1.1 to 1.5 according to barnea 1987
        gamma = 1.1

    # bubble rise velocity
    bubble_rise_vel = 1.53 * (grav * (rho_l - rho_g) * sigma / (rho_l ** 2)) ** (1 / 4)

    # lhs of condition
    angle_term = np.cos(beta) / (np.sin(beta) ** 2)
    # rhs of condition
    velocity_term = (
        ((3 / 4) * np.cos(np.pi / 4))
        * ((bubble_rise_vel ** 2) / grav)
        * ((lift_coef * gamma ** 2) / diam)
    )

    # if the angle is large enough to prevent bubbles from accumulating at the top
    return angle_term > velocity_term


def gas_void_fraction(u_gs, u_ls, liquid, gas, pipe, critical_value=0.25):
    """check if the gas void fraction is above the threshold
    at which bubbles turn into slug flow
        Taitel et al. 1980
    """
    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    grav = pipe.gravity
    beta = pipe.inclination

    # constant for convenience and readability
    k = 1.53 * np.sin(beta) * (grav * (rho_l - rho_g) * sigma / (rho_l ** 2)) ** (1 / 4)

    # terms for readability
    rhs_1 = u_gs * (1 - critical_value) / critical_value
    rhs_2 = k * (1 - critical_value)
    rhs = rhs_1 - rhs_2

    return u_ls > rhs
