""" calculate non dimensional numbers
"""
import numpy as np

from . import general


def reynolds(velocity, fluid, pipe):
    """calculate the reynolds number based on the fluid and its velocity
    in the pipe
    """
    return velocity * fluid.density * pipe.diameter / fluid.dynamic_viscosity


def lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe):
    """calculate the lockhardt martinelli number of the u_gs,
    u_ls combination
    """

    # split in terms for readability
    # numerator
    dpdx_ls = general.single_phase_dpdx(u_ls, liquid, pipe)

    # denominator
    dpdx_gs = general.single_phase_dpdx(u_gs, gas, pipe)

    lock_mart_number = np.sqrt(dpdx_ls / dpdx_gs)

    return lock_mart_number


def y_gravity(u_gs, u_ls, liquid, gas, pipe):
    """calculate the y value which represents the relative
    forces acting on the fluid in the flow direction due to
    gravity and pressure drop
    """
    # local variables
    rho_g = gas.density
    rho_l = liquid.density
    grav = pipe.gravity
    beta = pipe.inclination

    dpdx_gs = general.single_phase_dpdx(u_gs, gas, pipe)

    y_grav = (rho_l - rho_g) * grav * np.sin(beta) / dpdx_gs

    return y_grav
