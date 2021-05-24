"""
the functions that define the conditions for flow to be considered
annular flow
"""
import numpy as np
from scipy.optimize import newton
import general


def liquid_stability(u_gs, u_ls, liquid, gas, pipe):
    """find the locations at which the films stability condition is met"""

    # get the non dimensional values
    x_sqrd = general.lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe) ** 2
    y_grav = general.y_gravity(u_gs, u_ls, liquid, gas, pipe)

    # iterate to find the holdup
    initial_alpha_l = general.fluid_area_ratio(u_ls, liquid, pipe)
    liquid_holdup = newton(
        Calculate.liquid_instability, initial_alpha_l, args=(x_sqrd,)
    )

    # calculate the condition
    rhs = Calculate.equation16_barnea1987(liquid_holdup, x_sqrd)

    # the condition
    return y_grav < rhs


class Calculate:
    """
    namespace for the functions that are calculations and not maps
    """

    def liquid_instability(alpha_l, x_sqrd):
        """f(x)=0 formulation to find the equivalent annular liquid
        holdup at every u_gs, u_ls location
        barnea 1987 equation 15 and 16
        """
        # equation 15
        lhs = Calculate.equation15_barnea1987(alpha_l, x_sqrd)
        # equation 16
        rhs = Calculate.equation16_barnea1987(alpha_l, x_sqrd)
        return lhs - rhs

    def equation15_barnea1987(alpha_l, x_sqrd):
        """
        equation 15 in barnea 1987
        """
        term_1 = (1 + 75 * alpha_l) / (alpha_l * ((1 - alpha_l) ** (5 / 2)))
        term_2 = (1 / (alpha_l ** 3)) * x_sqrd
        y_grav = term_1 - term_2
        return y_grav

    def equation16_barnea1987(alpha_l, x_sqrd):
        """
        equation 16 in barnea 1987
        """
        numerator = 2 - (3 / 2) * alpha_l
        denominator = (alpha_l ** 3) * (1 - (3 / 2) * alpha_l)
        y_grav = (numerator / denominator) * x_sqrd
        return y_grav
