"""
define the constants of the two phases. Soft dependency on coolprop for some functionality
"""


class Gas:
    """
    define a class with gas properties
    """

    def __init__(self, density):
        self.density = density


class Liquid:
    """
    define a class with liquid properties
    """

    def __init__(self, density, bubble_surface_tension):
        self.density = density
        self.bubble_surface_tension = bubble_surface_tension


class Pipe:
    """
    define a class with pipe constants
    """

    def __init__(self, pipe_diameter, gravity=9.8):
        self.diameter = pipe_diameter
        self.gravity = gravity
