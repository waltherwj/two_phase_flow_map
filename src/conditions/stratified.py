"""
the functions that define the conditions for flow to be considered
stratified flow
"""
import numpy as np
from scipy.optimize import newton
import general


def wave_growth(u_gs, u_ls, liquid, gas, pipe):
    """calculate the conditon at which waves will not yet grow
    Taitel&Duckler 1976
    """

    # local variables for readability
    diam = pipe.diameter
    area = pipe.area

    # get the modified froude number
    froude = Calculate.modified_froude(u_gs, liquid, gas, pipe)

    # calculate the height corresponding to the area
    height = Calculate.sagitta_height(u_ls, liquid, pipe)

    # get the non-dimensional values
    # height ratio
    height_tilde = height / diam

    # areas
    area_ratio_gas = general.fluid_area_ratio(u_gs, gas, pipe)
    area_gas = area_ratio_gas * area
    area_gas_tilde = area_gas / (diam ** 2)

    # derivative
    r = diam / 2
    h = height
    dadh_l_tilde = (
        (2 * r ** 2) / np.sqrt(1 - (1 - 2 * h) ** 2)
        + 2 * r * np.sqrt(r ** 2 - (r - 2 * h * r) ** 2)
        - (2 * r * (r - 2 * h * r) ** 2) / np.sqrt(r ** 2 - (r - 2 * h * r) ** 2)
    )

    # velocity
    u_g_tilde = area / area_gas

    # divide lhs in terms for readability
    lhs_1 = 1 / ((1 - height_tilde) ** 2)
    lhs_2 = (u_g_tilde ** 2) * dadh_l_tilde / area_gas_tilde
    lhs = (froude ** 2) * lhs_1 * lhs_2

    # if lhs<1, waves are not yet forming
    return lhs < 1


class Calculate:
    """
    namespace for the methods that perform calculations as opposed to
    returning condition maps
    """

    @staticmethod
    def sagitta_height(velocity, fluid, pipe):
        """
        takes single fluid velocities and returns the equivalent
        height ratios
        """

        radius = pipe.diameter / 2

        # calculate the area ratio the fluid occupies
        area_ratio = general.fluid_area_ratio(velocity, fluid, pipe)
        valid_mask = area_ratio < 1 & area_ratio > 0

        # calculate the absolute area that corresponds to this
        area_fluid = area_ratio * pipe.area

        # iterate to find the height ratio that corresponds to this
        # the function
        area_function = (
            lambda h, r, area: area
            - (r ** 2) * np.arccos(1 - h / r)
            - (r - h) * np.sqrt(r ** 2 - (r - h) ** 2)
        )

        # initialize the array
        height = np.zeros_like(area_ratio)

        # use newton to find it
        initial_guess = np.ones_like(area_ratio[valid_mask]) * 0.25
        height[valid_mask] = newton(
            area_function, initial_guess, args=(radius, area_fluid[valid_mask])
        )

        # update the values that haven't been solved
        height[area_ratio >= 1] = pipe.diameter
        height[area_ratio <= 0] = 0

        return height

    @staticmethod
    def modified_froude(u_gs, liquid, gas, pipe):
        """calculate the modified froude modified by density ratio"""

        # local variables for readability
        rho_l = liquid.density
        rho_g = gas.density
        diam = pipe.diameter
        grav = pipe.gravity
        beta = pipe.inclination

        # terms for readability
        froude_1 = np.sqrt(rho_l / (rho_l - rho_g))
        froude_2 = u_gs / np.sqrt(diam * grav * np.cos(beta))
        froude = froude_1 * froude_2

        return froude
