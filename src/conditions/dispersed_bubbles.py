"""
the functions that define the conditions for flow to be considered
dispersed bubble flow
"""

import numpy as np
import general.friction_factor as friction_factor
import general
from fluids import Mix


# def gas_void_fraction(u_gs, u_ls, liquid, gas, pipe, critical_value=0.52):
#     """check if the gas void fraction is above the threshold
#     at which bubbles turn into slug flow
#         Taitel et al. 1980
#     """
#     # local variables for readability
#     rho_l = liquid.density
#     rho_g = gas.density
#     sigma = liquid.bubble_surface_tension
#     grav = pipe.gravity
#     beta = pipe.inclination

#     # constant for convenience and readability
#     k = 1.53 * np.sin(beta) * (grav * (rho_l - rho_g) * sigma / (rho_l ** 2)) ** (1 / 4)

#     # calculate the gas void fraction. Several results depending on the scenario
#     if k != 0:
#         mutable_term = np.sqrt(-4 * u_gs * k + (u_gs + k + u_ls) ** 2)
#         constant_term = u_gs + k + u_ls

#         # calculate the values
#         gas_void_frac = (constant_term - mutable_term) / (2 * k)

#         # check in which locations solutions are not valid, and replace with other
#         invalid_mask = (gas_void_frac < 0) & (gas_void_frac > 1)
#         gas_void_frac[invalid_mask] = (
#             constant_term[invalid_mask] + mutable_term[invalid_mask]
#         ) / (2 * k)
#     else:
#         gas_void_frac = u_gs / (u_ls + u_gs)

#     return gas_void_frac < critical_value


def gas_void_fraction(u_gs, u_ls, liquid, gas, pipe, critical_value=0.52):
    return general.fluid_area_ratio(u_gs, gas, pipe) < 0.52


def bubble_coalescence(u_gs, u_ls, liquid, gas, pipe):
    """check if the bubbles are large enough to coalesce into
    slugs with consideration for angle
     Barnea 1980
    """

    # local variables for readability
    roughness = pipe.roughness
    sigma = liquid.bubble_surface_tension
    rho_l = liquid.density
    diam = pipe.diameter

    # calculate the lhs
    deformed_bubble_size = Calculate.deformed_bubble_critical_size(liquid, gas, pipe)
    migration_to_top_size = Calculate.migration_to_top_critical_size(
        u_gs, u_ls, liquid, gas, pipe
    )
    # get a mask where the deformed bubble is the critical size
    mask = deformed_bubble_size < migration_to_top_size
    # make the array
    bubble_crit_diam = migration_to_top_size
    bubble_crit_diam[mask] = deformed_bubble_size

    # calculate the rhs

    # get the mixture values
    mix = Mix(u_gs, u_ls, liquid, gas, pipe)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)

    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)

    # mixture friction factor
    fric_mix = friction_factor.laminar_and_fang(reynolds_mix, roughness)

    # get the terms for readability
    rhs_1 = 0.725 + 4.15 * np.sqrt(u_gs / u_mix)
    rhs_2 = (sigma / rho_l) ** (3 / 5)
    rhs_3 = ((2 * fric_mix / diam) * (u_mix ** 3)) ** (-2 / 3)
    rhs = rhs_1 * rhs_2 * rhs_3

    # check if the bubbles are few/small enough that they won't coalesce
    return bubble_crit_diam < rhs


class Calculate:
    """
    namespace for the methods that perform calculations as opposed to
    returning condition maps
    """

    @staticmethod
    def deformed_bubble_critical_size(liquid, gas, pipe):
        """calculate the critical size of a bubble for it to become deformed
        eq. 6 Barnea 1987
        """
        # local variables for readability
        rho_l = liquid.density
        rho_g = gas.density
        sigma = liquid.bubble_surface_tension
        grav = pipe.gravity

        # critical bubble size
        diam_crit_deformed = 2 * np.sqrt(0.4 * sigma / ((rho_l - rho_g) * grav))

        return diam_crit_deformed

    @staticmethod
    def migration_to_top_critical_size(u_gs, u_ls, liquid, gas, pipe):
        """calculate the critical size below which a bubble can't travel
        to the upper part of pipe
        eq. 6 Barnea 1987
        """
        # get the mixture
        mix = Mix(u_gs, u_ls, liquid, gas, pipe)

        # local variables for readability
        rho_l = liquid.density
        rho_g = gas.density
        roughness = pipe.roughness
        grav = pipe.gravity
        beta = pipe.inclination

        # get the mixture velocity
        u_mix = mix.mixture_velocity(u_gs, u_ls)

        # mixture reynolds number
        reynolds_mix = general.reynolds(u_mix, mix, pipe)

        # mixture friction factor
        friction_mix = friction_factor.laminar_and_fang(reynolds_mix, roughness)

        # critical bubble size
        diam_crit_migration = (
            (3 / 8) * (rho_l / (rho_l - rho_g)) * (friction_mix * (u_mix ** 2))
        ) / (grav * np.cos(beta))

        return diam_crit_migration
