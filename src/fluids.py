"""
define the constants of the two phases. Soft dependency on coolprop for some functionality
"""
import numpy as np
from general import fluid_area_ratio


class Gas:
    """
    define a class with gas properties
    """

    def __init__(self, mass_flowrate, density, dynamic_viscosity):
        self.mass_flowrate = mass_flowrate
        self.density = density
        self.dynamic_viscosity = dynamic_viscosity


class Liquid:
    """
    define a class with liquid properties
    """

    def __init__(
        self, mass_flowrate, density, dynamic_viscosity, bubble_surface_tension
    ):
        self.mass_flowrate = mass_flowrate
        self.density = density
        self.bubble_surface_tension = bubble_surface_tension
        self.dynamic_viscosity = dynamic_viscosity


class Pipe:
    """
    define a class with pipe constants
    """

    def __init__(self, diameter, inclination, roughness, gravity=9.81):
        self.diameter = diameter
        self.gravity = gravity
        self.roughness = roughness
        self.inclination = inclination * np.pi / 180
        self.area = (diameter ** 2) * np.pi / 4


class Mix:
    """
    define a class with 2-phase mixed properties
    """

    def __init__(self, u_gs, u_ls, liquid, gas, pipe, gas_holdup=None):

        self.liquid = liquid
        self.gas = gas
        self.pipe = pipe
        self.gas_holdup = gas_holdup
        self.mass_flowrate = liquid.mass_flowrate + gas.mass_flowrate

        # the values that need to be calculated
        self.density = self.mixture_density(u_gs, u_ls)
        self.dynamic_viscosity = self.mixture_viscosity()

    @staticmethod
    def mixture_velocity(u_gs, u_ls):
        """calculate the mixture velocity"""
        return u_gs + u_ls

    def mixture_density(self, u_gs, u_ls):
        """calculate the equivalent mixture density"""

        # local variables for readability
        rho_l = self.liquid.density
        rho_g = self.gas.density

        # get the area ratios corresponding to superficial velocity values
        if self.gas_holdup is None:
            gas_area_ratio = u_gs / self.mixture_velocity(u_gs, u_ls)
        else:
            gas_area_ratio = self.gas_holdup

        # handle unreasonable values
        gas_area_ratio[gas_area_ratio >= 1] = 1 - 1e-12
        gas_area_ratio[gas_area_ratio <= 0] = 1e-12
        gas_area_ratio[np.isnan(gas_area_ratio)] = (
            u_gs / self.mixture_velocity(u_gs, u_ls)
        )[np.isnan(gas_area_ratio)]

        # the liquid ratio
        liq_area_ratio = 1 - gas_area_ratio

        # get the equivalent mixture density
        rho_mix = gas_area_ratio * rho_g + liq_area_ratio * rho_l

        return rho_mix

    def mixture_viscosity(self, model="owen"):
        """calculate the equivalent mixture density
        can have different models by selecting a keyword
        """

        # this is the preferred model according to N. Z. Aung and T. Yuwono, 2012
        if model == "owen":
            viscosity = self.liquid.dynamic_viscosity

        return viscosity
