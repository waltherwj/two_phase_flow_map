"""
the functions that define the conditions for flow to be considered
slug flow
"""
import numpy as np


def gas_void_fraction(ugs, uls, liquid, gas, pipe, critical_value=0.25):
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
        mutable_term = np.sqrt(-4 * ugs * k + (ugs + k + uls) ** 2)
        constant_term = ugs + k + uls

        # calculate the values
        gas_void_frac = (constant_term - mutable_term) / (2 * k)

        # check in which locations solutions are not valid, and replace with other
        invalid_mask = (gas_void_frac < 0) & (gas_void_frac > 1)
        gas_void_frac[invalid_mask] = (
            constant_term[invalid_mask] - mutable_term[invalid_mask]
        ) / (2 * k)
    else:
        gas_void_frac = ugs / (uls + ugs)

    return gas_void_frac > critical_value
