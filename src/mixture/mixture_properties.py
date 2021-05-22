""" module to calculate properties of a mixture
"""


def mixture_velocity(u_gs, u_ls):
    """calculate the mixture velocity"""
    return u_gs + u_ls


def mixture_density(u_gs, u_ls, liquid, gas, pipe):
    """calculate the equivalent mixture density"""
    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density

    # get the area ratios corresponding to superficial velocity values
    gas_area_ratio = fluid_area_ratio(u_gs, gas, pipe)
    liq_area_ratio = fluid_area_ratio(u_ls, liquid, pipe)

    # get the equivalent mixture density
    rho_mix = gas_area_ratio * rho_g + liq_area_ratio * rho_l

    return rho_mix


def mixture_viscosity(u_gs, u_ls, liquid, gas, pipe):
    """calculate the equivalent mixture density"""
    # local variables for readability
    rho_l = liquid.density
    rho_g = gas.density

    # get the area ratios corresponding to superficial velocity values
    gas_area_ratio = fluid_area_ratio(u_gs, gas, pipe)
    liq_area_ratio = fluid_area_ratio(u_ls, liquid, pipe)

    # get the equivalent mixture density
    rho_mix = gas_area_ratio * rho_g + liq_area_ratio * rho_l

    return rho_mix
