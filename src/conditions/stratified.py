"""
the functions that define the conditions for flow to be considered
stratified flow
"""
import numpy as np
from scipy.optimize import newton
import general
from general.general import Geometry
from general import friction_factor, non_dimensional


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
    return lhs - 1


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
        """calculate the equilibrium level of liquid at every
        single u_gs, u_ls pair
        """
        h_init = np.ones_like(u_gs) * 0.25
        h_equilibrium = newton(
            Calculate.equilibrium_equation, h_init, args=(u_gs, u_ls, liquid, gas, pipe)
        )

        return h_equilibrium

    @staticmethod
    def equilibrium_equation(height_tilde, u_gs, u_ls, liquid, gas, pipe):
        """f(x) = 0 formulation to iterate and find the equilibrium level of liquid at every
        single u_gs, u_ls pair
        """

        # get the geometry values for these heights
        geom = Geometry(height_tilde, absolute=True, pipe=pipe)

        # local variables for readability and correspondence to the equation
        rho_l = liquid.density
        rho_g = gas.density
        mu_l = liquid.dynamic_viscosity
        mu_g = gas.dynamic_viscosity
        roughness = pipe.roughness
        beta = pipe.inclination
        grav = pipe.gravity
        s_gas = geom.perim_gas
        s_liq = geom.perim_liq
        s_interf = geom.perim_interf
        a_gas = geom.a_gas
        a_liq = geom.a_liq

        # average velocities
        u_g_mean = u_gs / a_gas
        u_l_mean = u_ls / a_liq
        u_i_mean = u_l_mean

        # hydraulic diameters
        hydr_diam_liq = 4 * a_liq / s_liq
        hydr_diam_gas = 4 * a_gas / (s_liq + s_interf)

        # reynold numbers
        reynolds_l = u_l_mean * rho_l * hydr_diam_liq / mu_l
        reynolds_g = u_g_mean * rho_g * hydr_diam_gas / mu_g

        # friction factors
        fric_liq = friction_factor.laminar_and_fang(reynolds_l, roughness)
        fric_gas = friction_factor.laminar_and_fang(reynolds_g, roughness)
        fric_interf = fric_gas

        # shear stresses
        shear_gas = fric_gas * rho_g * (u_g_mean ** 2) / 2
        shear_liq = fric_liq * rho_l * (u_l_mean ** 2) / 2
        shear_interf = fric_interf * rho_g * ((u_g_mean - u_i_mean) ** 2) / 2

        # balance equation
        gas_term = shear_gas * s_gas / a_gas
        liq_term = shear_liq * s_liq / a_liq
        interf_term = shear_interf * s_interf * ((1 / a_gas) + (1 / a_liq))
        grav_term = (rho_l - rho_g) * grav * np.sin(beta)

        return gas_term - liq_term + interf_term + grav_term

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
