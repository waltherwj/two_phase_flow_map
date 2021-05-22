"""
the functions that define the conditions for flow to be considered
bubbly flow
"""


def taylor_bubble_exceeds(liquid, gas, pipe):
    """check if taylor bubble velocity exceeds bubble velocity
    which is necessary for bubbly flow to exist
        Taitel et al. 1980
    """
    rho_l = liquid.density
    rho_g = gas.density
    sigma = liquid.bubble_surface_tension
    gravity = pipe.gravity

    rhs = 19 * (rho_l - rho_g)
