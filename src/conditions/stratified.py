"""
the functions that define the conditions for flow to be considered
stratified flow
"""
import numpy as np
from scipy.optimize import newton
import general
from general.general import Geometry
from general import friction_factor, non_dimensional


def wave_growth(u_gs, u_ls, liquid, gas, pipe):
    """check the conditon at which waves will not yet grow
    Taitel&Duckler 1976
    """

    # local variables for readability
    diam = pipe.diameter
    area = pipe.area

    # get the modified froude number
    froude = Calculate.modified_froude(u_gs, liquid, gas, pipe)

    # height
    height = Calculate.sagitta_absolute_height(u_ls, liquid, pipe)

    # get the non-dimensional values
    # height ratio
    height_tilde = height / diam

    # non dimensional values
    tilde = Geometry(height_tilde)

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
    return lhs < 1


class Calculate:
    """
    namespace for the methods that perform calculations as opposed to
    returning condition maps
    """

    def sagitta_absolute_height(velocity, fluid, pipe):
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
        """calculate the equilibrium level of liquid at every
        single u_gs, u_ls pair
        """
        h_init = np.ones_like(u_gs) * 0.25
        h_equilibrium = newton(
            Calculate.equilibrium_equation, h_init, args=(u_gs, u_ls, liquid, gas, pipe)
        )

        return h_equilibrium

    @staticmethod
    def equilibrium_equation(u_gs, u_ls, liquid, gas, pipe):
        """f(x) = 0 formulation to iterate and find the equilibrium level of liquid at every
        single u_gs, u_ls pair
        """

        # local variables for readability and correspondence to the equation
        rho_l = liquid.density
        rho_g = gas.density
        mu_l = liquid.dynamic_viscosity
        mu_g = gas.dynamic_viscosity

        # get the geometry values for these heights
        height_tilde = Calculate.sagitta_absolute_height(u_ls, liquid, pipe)
        tilde = Geometry(height_tilde, non_dimensional=True, pipe=pipe)

        # get the non dimensional numbers
        x_sqrd = non_dimensional.lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe) ** 2
        y_grav = non_dimensional.y_gravity(u_gs, u_ls, liquid, gas, pipe)

        # get the single fluid reynolds numbers from the non dimensional values
        reynolds_ls = (
            (rho_l * u_ls * pipe.diameter / mu_l) * tilde.vel_l * tilde.hydr_diam_l
        )
        reynolds_gs = (
            (rho_g * u_gs * pipe.diameter / mu_g) * tilde.vel_g * tilde.hydr_diam_g
        )

        # get the exponents for the blasius friction factor
        m = Calculate.blasius_exponents(reynolds_gs)
        n = Calculate.blasius_exponents(reynolds_ls)

        # balance equation
        gas_term = (
            ((tilde.hydr_diam_g * tilde.vel_g) ** (-m))
            * (tilde.vel_g ** 2)
            * (
                (tilde.perim_g / tilde.area_g)
                + (tilde.perim_interf / tilde.area_l)
                + (tilde.perim_interf / tilde.area_g)
            )
        )
        liq_term = (
            x_sqrd
            * ((tilde.vel_l * tilde.hydr_diam_l) ** (-n))
            * (tilde.vel_l ** 2)
            * (tilde.perim_l / tilde.area_l)
        )
        grav_term = 4 * y_grav

        return gas_term - liq_term - grav_term > 0

    @staticmethod
    def blasius_exponents(reynolds):
        """calculates the exponents related to the blasius equation"""
        exponents = np.empty_like(reynolds)
        exponents[reynolds > 2300] = 0.2
        exponents[reynolds <= 2300] = 1

        return exponents

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
