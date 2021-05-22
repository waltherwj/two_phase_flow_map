"""functions which are too broad to fit in any other category
and are used by several different other modules
"""


from conditions.dispersed_bubbles import gas_void_fraction


def fluid_area_ratio(velocity, fluid, pipe):
    """calculate the the fluid are ratio based on the fluid velocity
    and properties.

    Can be used to calculate gas void fraction
    """
    area_gas = fluid.mass_flux / (fluid.density * velocity)

    area_ratio = area_gas / pipe.area

    return area_ratio


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
