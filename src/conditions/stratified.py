"""
the functions that define the conditions for flow to be considered
stratified flow
"""
import numpy as np
from scipy.optimize import newton
import general


def wave_growth(u_gs, u_ls, height, liquid, gas, pipe):
    """check the conditon at which waves will not yet grow
    Taitel&Duckler 1976
    """

    # local variables for readability
    diam = pipe.diameter
    area = pipe.area

    # get the modified froude number
    froude = Calculate.modified_froude(u_gs, liquid, gas, pipe)

    # get the non-dimensional values
    # height ratio
    height_tilde = height / diam

    # non dimensional values
    tilde = NonDimGeom(height_tilde)

    # areas
    area_ratio_gas = general.fluid_area_ratio(u_gs, gas, pipe)
    area_gas = area_ratio_gas * area
    area_gas_tilde = area_gas / (diam ** 2)

    # derivative
    dadh_l_tilde = ((2 - 2 * height_tilde) * height_tilde) / np.sqrt(
        (1 - height_tilde) * height_tilde
    )

    # velocity
    u_g_tilde = area / area_gas

    # divide lhs in terms for readability
    lhs_1 = 1 / ((1 - height_tilde) ** 2)
    lhs_2 = (u_g_tilde ** 2) * dadh_l_tilde / area_gas_tilde
    lhs = (froude ** 2) * lhs_1 * lhs_2

    # if lhs<1, waves are not yet forming
    return lhs - 1


class NonDimGeom:
    """calculates and stores the non dimensional geometric values related to the
    liquid height
    """

    def __init__(self, height_ratio):

        # height
        height = self.sagitta_absolute_height(u_ls, liquid, pipe)
        self.height = height / pipe.diameter

        # areas
        area_gas_ratio = general.fluid_area_ratio(u_gs, gas, pipe)
        area_liq_ratio = general.fluid_area_ratio(u_ls, liquid, pipe)
        area_gas = area_gas_ratio * pipe.area
        area_liq = area_liq_ratio * pipe.area
        self.area_gas = area_gas / (pipe.diameter ** 2)
        self.area_liq = area_liq / (pipe.diameter ** 2)


class Calculate:
    """
    namespace for the methods that perform calculations as opposed to
    returning condition maps
    """

    def sagitta_absolute_height(self, velocity, fluid, pipe):
        """
        takes single fluid velocities and returns the equivalent
        height ratios
        """

        radius = pipe.diameter / 2

        # calculate the area ratio the fluid occupies
        area_ratio = general.fluid_area_ratio(velocity, fluid, pipe)
        valid_mask = (area_ratio < 1) & (area_ratio > 0)

        # calculate the absolute area that corresponds to this
        area_fluid = area_ratio * pipe.area

        # iterate to find the angle theta that corresponds to this
        # the function
        area_function = lambda theta, r, area: area - ((r ** 2) / 2) * (
            theta - np.sin(theta)
        )

        # initialize the array
        theta = np.zeros_like(area_ratio)

        # use newton to find it
        initial_guess = np.ones_like(area_ratio[valid_mask]) * np.pi / 4
        theta[valid_mask] = newton(
            area_function, initial_guess, args=(radius, area_fluid[valid_mask])
        )

        # update the values that haven't been solved
        theta[area_ratio >= 1] = 2 * np.pi
        theta[area_ratio <= 0] = 0

        # find the heights from the theta
        height = radius * (1 - np.cos(theta / 2))

        return height

    @staticmethod
    def equilibrium_level(u_gs, u_ls, liquid, gas, pipe):
        """
        calculate the equilibrium level of liquid at every
        single u_gs, u_ls pair
        """

        pass

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
