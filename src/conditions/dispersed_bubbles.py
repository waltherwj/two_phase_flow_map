"""
the functions that define the conditions for flow to be considered
dispersed bubble flow
"""

from general.non_dimensional import reynolds
import numpy as np
import general.friction_factor as friction_factor
import general
from fluids import Mix


def gas_void_fraction(u_gs, u_ls, critical_value=0.35):
    """condition for gas void fraction
    the realistic gas void fraction is not 0.52,
    as higher gas void fractions than 0.35 are rare
    TODO find citation
    """
    return u_ls > u_gs * (1 - critical_value) / critical_value


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
    # make the array
    bubble_crit_diam = np.minimum(migration_to_top_size, deformed_bubble_size)

    # get the mixture values
    mix = Mix(u_gs, u_ls, liquid, gas, pipe)
    # get the mixture velocity
    u_mix = mix.mixture_velocity(u_gs, u_ls)

    # mixture reynolds number
    reynolds_mix = general.reynolds(u_mix, mix, pipe)

    # mixture friction factor
    fric_mix = friction_factor.fang(reynolds_mix, roughness)

    # get the terms for readability
    rhs_1 = 0.725 + 4.15 * np.sqrt(u_gs / u_mix)
    rhs_2 = (sigma / rho_l) ** (3 / 5)
    rhs_3 = ((2 * fric_mix / diam) * (u_mix ** 3)) ** (-2 / 5)
    rhs = rhs_1 * rhs_2 * rhs_3

    # check if the bubbles are few/small enough that they won't coalesce
    return bubble_crit_diam > rhs


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
        friction_mix = friction_factor.fang(reynolds_mix, roughness)

        # critical bubble size
        diam_crit_migration = (
            (3 / 8) * (rho_l / (rho_l - rho_g)) * (friction_mix * (u_mix ** 2))
        ) / (grav * np.cos(beta))

        return diam_crit_migration
