""" calculate non dimensional numbers
"""


def reynolds(velocity, fluid, pipe):
    """calculate the reynolds number based on the fluid and its velocity
    in the pipe
    """
    return velocity * fluid.density * pipe.diameter / fluid.dynamic_viscosity
