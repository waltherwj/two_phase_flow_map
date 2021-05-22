"""
the functions that define the conditions for flow to be considered
slug flow
"""
import numpy as np


def gas_void_fraction(u_gs, u_ls, liquid, gas, pipe, critical_value=0.25):
    """check if the gas void fraction is above the threshold
    at which bubbles turn into slug flow
        Taitel et al. 1980
    """
    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    grav = pipe.gravity
    beta = pipe.inclination

    # constant for convenience and readability
    k = 1.53 * np.sin(beta) * (grav * (rho_l - rho_g) * sigma / (rho_l ** 2)) ** (1 / 4)

    # calculate the gas void fraction. Several results depending on the scenario
    if k != 0:
        mutable_term = np.sqrt(-4 * u_gs * k + (u_gs + k + u_ls) ** 2)
        constant_term = u_gs + k + u_ls

        # calculate the values
        gas_void_frac = (constant_term - mutable_term) / (2 * k)

        # check in which locations solutions are not valid, and replace with other
        invalid_mask = (gas_void_frac < 0) & (gas_void_frac > 1)
        gas_void_frac[invalid_mask] = (
            constant_term[invalid_mask] + mutable_term[invalid_mask]
        ) / (2 * k)
    else:
        gas_void_frac = u_gs / (u_ls + u_gs)

    return gas_void_frac > critical_value
