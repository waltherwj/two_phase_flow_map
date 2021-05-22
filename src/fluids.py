"""
define the constants of the two phases. Soft dependency on coolprop for some functionality
"""
import numpy as np


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

    def __init__(self, diameter, inclination, gravity=9.81):
        self.diameter = diameter
        self.gravity = gravity
        self.inclination = inclination * np.pi / 180
