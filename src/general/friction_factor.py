""" Provides various functions for calculating the Darcy friction 
factor based on Colebrook-White approximations.

adapted from 
https://github.com/IMEConsultants/colebrook/blob/master/colebrook/colebrook.py
 to be compatible with numpy arrays. Under GNU 3 License
"""


import numpy as np


def sjFriction(reynolds, roughness):
    """
    Model: Swamee, Jain
    Year: 1976
    Paper: https://cedb.asce.org/CEDBsearch/record.jsp?dockey=0006693
    Suitable Range:
        5000 < Reyonolds < 10^8 and e/D = 0.00001 - 0.5
    """
    friction = 0.25 / (np.log10(roughness / 3.7 + 5.74 / (reynolds ** 0.9))) ** 2
    return friction


def bntFriction(reynolds, roughness):
    """
    Model: Bellos, Nalbantis, Tsakris
    Year: 2018
    Paper: https://www.sciencedirect.com/science/article/pii/S0029549311000173
    Suitable Range:
        All Flow Regimes
    """
    reynolds = np.array(reynolds)
    roughness = np.array(roughness)
    friction = np.empty_like(reynolds)

    # laminar + transition
    low_re_mask = reynolds < 3000
    friction[low_re_mask] = 0.316 / (reynolds[low_re_mask] ** (1 / 4))

    # transition + turbulent
    roughness_term = 0.234 * (roughness ** 1.1007)[~low_re_mask]
    reynolds_term_neg = -60.525 * (reynolds ** 1.1105)[~low_re_mask]
    reynolds_term_pos = 56.291 * (reynolds ** 1.0712)[~low_re_mask]
    friction[~low_re_mask] = 1.613 * (
        np.log(roughness_term + reynolds_term_neg + reynolds_term_pos) ** (-2)
    )

    return friction


def fngFriction(reynolds, roughness):
    """
    Model: Fang
    Year: 2011
    Paper: https://www.sciencedirect.com/science/article/pii/S0029549311000173
    Suitable Range:
        Reynolds > 2300 (I.E. Turbulent and Transition Range only)
    """

    friction = (
        1.613
        * (
            np.log(
                0.234 * roughness ** 1.1007
                - 60.525 / reynolds ** 1.1105
                + 56.291 / reynolds ** 1.0712
            )
        )
        ** -2
    )
    return friction


def eptFriction(reynolds, roughness):
    """
    Model: Evangelides, Papaevangelou, Tzimopoulos
    Year: 2010
    Paper: http://blogs.sch.gr/geopapaevan/files/2010/07/full-paper_pre1128act.pdf
    Suitable Range:
        Reynolds > 2300 (I.E. Turbulent and Transition Range only)
    Rel
    """
    friction = (0.2479 - 0.0000947 * (7 - np.log10(reynolds)) ** 4) / (
        np.log10(roughness / 3.615 + 7.366 / reynolds ** 0.9142)
    ) ** 2
    return friction


def akFriction(reynolds, roughness):
    """
    Model: Avci, Kargoz
    Year: 2009
    Paper: http://dx.doi.org/10.1115/1.3129132
    Suitable Range:
        Reynolds > 2300 (I.E. Turbulent and Transition Range only)
    """
    friction = (
        6.4
        / (
            np.log(reynolds)
            - np.log(1 + 0.01 * reynolds * roughness * (1 + 10 * np.sqrt(roughness)))
        )
        ** 2.4
    )
    return friction


def bkcFriction(reynolds, roughness):
    """
    Model: Brkic
    Year: 2011
    Paper: https://doi.org/10.1080%2F10916461003620453
    Suitable Range:
        Reynolds > 2300 (I.E. Turbulent and Transition Range only)
    """
    beta = np.log(
        reynolds / (1.816 * np.log((1.1 * reynolds) / np.log(1 + 1.1 * reynolds)))
    )
    friction = (-2 * np.log10((2.18 * beta) / reynolds + roughness / 3.71)) ** -2
    return friction


if __name__ == "__main__":
    reynolds = np.random.rand() * 5000
    roughness = np.random.rand() * 0.05
    table_format = "|{:<15} | {:<7}|"
    fric_headers = ["METHOD", "VALUE"]

    factors = []
    factors.append(["sjFriction", sjFriction(reynolds, roughness)])
    factors.append(["bntFriction", bntFriction(reynolds, roughness)])
    factors.append(["fngFriction", fngFriction(reynolds, roughness)])
    factors.append(["akFriction", akFriction(reynolds, roughness)])
    factors.append(["bkcFriction", bkcFriction(reynolds, roughness)])
    factors.append(["eptFriction", eptFriction(reynolds, roughness)])
    print(
        "Ensure values are within range of applicability for equations (specifically around transition and laminar region)!"
    )
    print(table_format.format(*fric_headers))
    print("---------------------------")
    for row in factors:
        print(table_format.format(*row))
    print("DISCLAIMER: Use secondary verification. No guarantee of accuracy")
