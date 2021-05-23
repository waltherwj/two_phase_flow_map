"""functions which are too broad to fit in any other category
and are used by several different other modules
"""
from . import non_dimensional, friction_factor
import numpy as np


def fluid_area_ratio(velocity, fluid, pipe):
    """calculate the the fluid are ratio based on the fluid velocity
    and properties.

    Can be used to calculate gas void fraction
    """
    area_gas = fluid.mass_flowrate / (fluid.density * velocity)

    area_ratio = area_gas / pipe.area

    return area_ratio


def single_phase_dpdx(velocity, fluid, pipe):
    """get the dpdx of one phase flowing alone in the pipe"""
    # local variables
    rho = fluid.density
    roughness = pipe.roughness
    diam = pipe.diameter

    # get reynolds
    reynolds = non_dimensional.reynolds(velocity, fluid, pipe)

    # get friction factor
    fric = friction_factor.fang(reynolds, roughness)

    dpdx_s = (4 / diam) * fric * rho * (velocity ** 2) / 2

    return dpdx_s


class Geometry:
    """calculates and stores the geometric values related to the
    liquid height ratio
    """

    def __init__(
        self, height_ratio, non_dimensional=False, pipe=None, u_gs=None, u_ls=None
    ):

        if not (non_dimensional):
            diam = pipe.diameter
        else:
            diam = 1

        self.height_ratio = height_ratio
        self.height_ratio[height_ratio > 1] = 1
        self.height_ratio[height_ratio < 0] = 0
        # dummy constant for readability
        self.var = 2 * self.height_ratio - 1

        # the non dimensional constants
        # areas
        self.area_l = self.area_liq() * diam ** 2
        self.area_g = self.area_gas() * diam ** 2
        # perimeters
        self.perim_l = self.perimeter_liq() * diam
        self.perim_g = self.perimeter_gas() * diam
        self.perim_interf = self.perimeter_interface() * diam
        # hydraulic diameters
        self.hydr_diam_l = self.hydraulic_diameter(fluid="liq")
        self.hydr_diam_g = self.hydraulic_diameter(fluid="gas")
        # velocities
        if non_dimensional:
            self.vel_l = self.velocity_liq()
            self.vel_g = self.velocity_gas()
        elif (u_gs is not None) and (u_ls is not None):
            self.vel_l = self.velocity_liq() * u_ls
            self.vel_g = self.velocity_gas() * u_gs
        else:
            raise ValueError(
                "Geometry is not nondimensionalized."
                + " Need u_gs and u_ls for velocity calculation"
            )

    def hydraulic_diameter(self, fluid="liq"):
        """
        calculate the hydraulic diameter of the fluid
        """
        if fluid == "liq":
            return 4 * self.area_l / self.perim_l
        elif fluid == "gas":
            return 4 * self.area_g / (self.perim_g + self.perim_interf)
        else:
            raise RuntimeError(
                "hydraulic diameter needs <fluid> argument to be one of ['liq', 'gas']"
            )

    def area_liq(self):
        """
        calculate area_tilde, a non dimensionalized version of the area
        """
        var = self.var

        # terms for readability
        term_1 = np.pi - np.arccos(var)
        term_2 = var * np.sqrt(1 - var ** 2)
        liq_area_tilde = 0.25 * (term_1 + term_2)
        return liq_area_tilde

    def area_gas(self):
        """
        calculate area_tilde, a non dimensionalized version of the area
        """
        var = self.var

        # terms for readability
        term_1 = np.arccos(var)
        term_2 = var * np.sqrt(1 - var ** 2)
        gas_area_tilde = 0.25 * (term_1 - term_2)
        return gas_area_tilde

    def perimeter_liq(self):
        """
        calculate perimeter_tilde, a non dimensionalized version
        of the perimeter
        """
        return np.pi - np.arccos(self.var)

    def perimeter_gas(self):
        """
        calculate perimeter_tilde, a non dimensionalized version
        of the perimeter
        """
        return np.arccos(self.var)

    def perimeter_interface(self):
        """
        calculate perimeter_tilde, a non dimensionalized version
        of the perimeter
        """
        return np.sqrt(1 - self.var ** 2)

    def velocity_liq(self):
        """
        calculate velocity_tilde, a non dimensionalized version
        of the velocity
        """
        return (np.pi / 4) / self.area_liq()

    def velocity_gas(self):
        """
        calculate velocity_tilde, a non dimensionalized version
        of the velocity
        """
        return (np.pi / 4) / self.area_gas()


def sagitta_absolute_height(velocity, fluid, pipe):
    """
    takes single fluid velocities and returns the equivalent
    height
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
