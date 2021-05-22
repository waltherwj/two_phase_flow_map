""" Provides various functions for calculating the Darcy friction 
factor based on Colebrook-White approximations.

adapted from 
https://github.com/IMEConsultants/colebrook/blob/master/colebrook/colebrook.py
 to be compatible with numpy arrays. Under GNU 3 License
"""


import numpy as np


def laminar_and_fang(reynolds, roughness):
    """combine fang 2011 model for transition and turbulent with the laminar model"""

    friction = np.empty_like(reynolds)
    reynolds = np.array(reynolds)

    turbulent = reynolds >= 2300
    # laminar
    friction[~turbulent] = laminar(reynolds[~turbulent])

    # turbulent
    friction[turbulent] = fang(reynolds[turbulent], roughness)

    return friction


def laminar(reynolds):
    """calculate laminar friction factor"""
    return reynolds / 64


def turb_swamee(reynolds, roughness):
    """
    Model: Swamee, Jain
    Year: 1976
    Paper: https://cedb.asce.org/CEDBsearch/record.jsp?dockey=0006693
    Suitable Range:
        5000 < Reyonolds < 10^8 and e/D = 0.00001 - 0.5
    """
    friction = 0.25 / (np.log10(roughness / 3.7 + 5.74 / (reynolds ** 0.9))) ** 2
    return friction


def bnt(reynolds, roughness):
    """
    Model: Bellos, Nalbantis, Tsakris
    Year: 2018
    Paper
    Suitable Range:
        All Flow Regimes, FREE SURFACE FLOW
    """
    inv_roughness = 1 / roughness
    param_a = 1 / (1 + (reynolds / 2712) ** 8.4)
    param_b = 1 / (1 + (reynolds / (150 * inv_roughness)) ** 1.8)
    exponent_a = 2 * (param_a - 1) * param_b
    exponent_b = 2 * (param_a - 1) * (1 - param_b)
    friction = (
        (64 / reynolds) ** param_a
        * (0.75 * np.log(reynolds / 5.37)) ** exponent_a
        * (0.88 * np.log(6.82 * inv_roughness)) ** exponent_b
    )
    return friction


def fang(reynolds, roughness):
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


def ept(reynolds, roughness):
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


def ak(reynolds, roughness):
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


def bkc(reynolds, roughness):
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
    print(reynolds, roughness)
    table_format = "|{:<15} | {:<7}|"
    fric_headers = ["METHOD", "VALUE"]

    factors = []
    factors.append(["turb_swamee", turb_swamee(reynolds, roughness)])
    factors.append(["bnt", bnt(reynolds, roughness)])
    factors.append(["fng", fang(reynolds, roughness)])
    factors.append(["cfng", laminar_and_fang(reynolds, roughness)])
    factors.append(["ak", ak(reynolds, roughness)])
    factors.append(["bkc", bkc(reynolds, roughness)])
    factors.append(["ept", ept(reynolds, roughness)])
    print(
        "Ensure values are within range of applicability for equations (specifically around transition and laminar region)!"
    )
    print(table_format.format(*fric_headers))
    print("---------------------------")
    for row in factors:
        print(table_format.format(*row))
    print("DISCLAIMER: Use secondary verification. No guarantee of accuracy")