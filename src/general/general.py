"""functions which are too broad to fit in any other category
and are used by several different other modules
"""


def fluid_area_ratio(velocity, fluid, pipe):
    """calculate the the fluid are ratio based on the fluid velocity
    and properties.

    Can be used to calculate gas void fraction
    """
    area_gas = fluid.mass_flowrate / (fluid.density * velocity)

    area_ratio = area_gas / pipe.area

    return area_ratio
