"""
the functions that define the conditions for flow to be considered
bubbly flow
"""


def taylor_bubble_exceeds(constants):
    """check if taylor bubble velocity exceeds bubble velocity
    which is necessary for bubbly flow to exist
        Taitel et al. 1980
    """
    liq_dens = constants.liq_dens

    rhs = 19 * (liq_dens - gas_dens)
