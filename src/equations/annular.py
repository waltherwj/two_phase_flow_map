"""
the functions that to calculate the conditions for flow to be considered
annular flow
"""


def liquid_instability(alpha_l, y_grav, x_sqrd):
    """f(x)=0 formulation to find the equivalent annular liquid
    holdup at every u_gs, u_ls location
    barnea 1987 equation 15
    """
    # equation 15
    rhs = equation15_barnea1987(alpha_l, x_sqrd)

    return y_grav - rhs


def equation15_barnea1987(alpha_l, x_sqrd):
    """
    equation 15 in barnea 1987
    """
    term_1 = (1 + 75 * alpha_l) / (alpha_l * ((1 - alpha_l) ** (5 / 2)))
    term_2 = (1 / (alpha_l ** 3)) * x_sqrd
    y_grav = term_1 - term_2
    return y_grav


def equation16_barnea1987(alpha_l, x_sqrd):
    """
    equation 16 in barnea 1987
    """
    numerator = 2 - (3 / 2) * alpha_l
    denominator = (alpha_l ** 3) * (1 - (3 / 2) * alpha_l)
    y_grav = (numerator / denominator) * x_sqrd
    return y_grav
