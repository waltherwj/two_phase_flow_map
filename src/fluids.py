"""
define the fluid constants. Soft dependency on coolprop for some functionality
"""


class Fluids:
    """
    define a class fluid properties as an easy way to get values
    for the constants
    """

    def __init__(
        self,
        pipe_diameter,
        liq_density,
        gas_density,
        bubble_surface_tension,
        gravity=9.8,
    ):

        self.pipe_diameter = pipe_diameter
        self.liq_density = liq_density
        self.gas_density = gas_density
        self.bubble_surface_tension = bubble_surface_tension
        self.gravity = gravity
