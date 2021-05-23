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
    fric = friction_factor.laminar_and_fang(reynolds, roughness)

    dpdx_s = (4 / diam) * fric * rho * (velocity ** 2) / 2

    return dpdx_s


class Geometry:
    """calculates and stores the geometric values related to the
    liquid height ratio
    """

    def __init__(self, height_ratio, absolute=False, pipe=None):

        if absolute:
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
        self.a_liq = self.area_liq() * diam ** 2
        self.a_gas = self.area_gas() * diam ** 2
        # perimeters
        self.perim_liq = self.perimeter_liq() * diam
        self.perim_gas = self.perimeter_gas() * diam
        self.perim_interf = self.perimeter_interface() * diam
        if not (absolute):
            # velocities
            self.vel_liq = self.velocity_liq()
            self.vel_gas = self.velocity_gas()
        # hydraulic diameters

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
        term_1 = np.pi - np.arccos(var)
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
