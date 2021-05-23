"""functions which are too broad to fit in any other category
and are used by several different other modules
"""
from . import non_dimensional, friction_factor


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
