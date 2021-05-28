"""
the functions that define the conditions for flow to be considered
annular flow
"""
import numpy as np
from scipy.optimize import newton
import general
import equations


def liquid_stability(u_gs, u_ls, liquid, gas, pipe):
    """find the locations at which the films stability condition is met"""

    # get the non dimensional values
    x_sqrd = general.lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe) ** 2
    y_grav = general.y_gravity(u_gs, u_ls, liquid, gas, pipe)

    # iterate to find the holdup
    initial_alpha_l = 1 - u_gs / (u_ls + u_gs)
    liquid_holdup = newton(
        equations.annular.liquid_instability,
        initial_alpha_l,
        args=(
            y_grav,
            x_sqrd,
        ),
    )

    # get rid of nonsensical values
    liquid_holdup[(liquid_holdup < 0) | (liquid_holdup > 1)]
    # calculate the condition
    rhs = equations.annular.equation16_barnea1987(liquid_holdup, x_sqrd)

    # the condition
    return y_grav < rhs


def gas_core_blockage(u_gs, u_ls, liquid, gas, pipe):
    """
    find the locations in which the gas core is expercted to not be blocked
    """
    # get the non dimensional values
    x_sqrd = general.lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe) ** 2
    y_grav = general.y_gravity(u_gs, u_ls, liquid, gas, pipe)

    # iterate to find the holdup
    initial_alpha_l = 1 - u_gs / (u_ls + u_gs)
    liquid_holdup = newton(
        equations.annular.liquid_instability,
        initial_alpha_l,
        args=(
            y_grav,
            x_sqrd,
        ),
    )
    # get rid of nonsensical values
    liquid_holdup[(liquid_holdup < 0) | (liquid_holdup > 1)] = np.nan
    r_sm = 0.48

    return liquid_holdup / r_sm < 0.5
