"""
the functions that define the conditions for flow to be considered
stratified flow
"""
import numpy as np
from scipy.optimize import newton
from general.general import Geometry
from general import friction_factor, non_dimensional
import equations


def equilibrium_equation(u_gs, u_ls, liquid, gas, pipe):
    """find if it is stratified or not at the critical height
        based on the the equilibrium level of liquid at every
    single u_gs, u_ls pair
    """

    # local variables for readability and correspondence to the equation
    rho_l = liquid.density
    rho_g = gas.density
    mu_l = liquid.dynamic_viscosity
    mu_g = gas.dynamic_viscosity
    roughness = pipe.roughness

    # Set inital guess that solves from the top
    # the solution space is difficult
    # TODO find how to set a good initial guess
    height_initial = np.ones_like(u_gs) * 0.95

    # the critical heights at which waves would start to grow
    height_tilde = newton(
        equations.stratified.wave_growth, height_initial, args=(u_gs, liquid, gas, pipe)
    )

    # non dimensional
    tilde = Geometry(height_tilde, non_dimensional=True)
    # dimensional
    geom = Geometry(
        height_tilde, non_dimensional=False, pipe=pipe, u_gs=u_gs, u_ls=u_ls
    )

    # get the non dimensional numbers
    x_sqrd = non_dimensional.lockhart_martinelli(u_gs, u_ls, liquid, gas, pipe) ** 2
    y_grav = non_dimensional.y_gravity(u_gs, u_ls, liquid, gas, pipe)

    # get the single fluid reynolds numbers from the non dimensional values
    reynolds_ls = rho_l * u_ls * pipe.diameter / mu_l
    reynolds_gs = rho_g * u_gs * pipe.diameter / mu_g

    # get the actual fluid average reynolds
    reynolds_l_actual = rho_l * geom.vel_l * geom.hydr_diam_l / mu_l
    reynolds_g_actual = rho_l * geom.vel_g * geom.hydr_diam_g / mu_l

    # get the friction factors for all the fluids
    # actual
    friction_g = friction_factor.niazkar_and_churchill(reynolds_g_actual, roughness)
    friction_l = friction_factor.niazkar_and_churchill(reynolds_l_actual, roughness)
    # single fluid
    friction_gs = friction_factor.niazkar_and_churchill(reynolds_gs, roughness)
    friction_ls = friction_factor.niazkar_and_churchill(reynolds_ls, roughness)

    # balance equation adapted from dimensional version of Taitel 1976
    gas_term = (
        (friction_g / friction_gs)
        * (tilde.vel_g ** 2)
        * (
            (tilde.perim_g / tilde.area_g)
            + (tilde.perim_interf / tilde.area_l)
            + (tilde.perim_interf / tilde.area_g)
        )
    )
    liq_term = (
        x_sqrd
        * (friction_l / friction_ls)
        * (tilde.vel_l ** 2)
        * (tilde.perim_l / tilde.area_l)
    )
    grav_term = 4 * y_grav

    # this indicates the area where waves start to grow and the
    # equilibrium equation is not met
    return gas_term - liq_term - grav_term > 0


def too_steep_for_stratified(u_gs, u_ls, liquid, gas, pipe):
    """
    check if the liquid velocity is so high that in steep inclination it
    ends up tearing droplets apart from the wavy turbulent interface resulting
    in annular flow
    Barnea 1987 eq 13-14
    """
    rho_l = liquid.density
    mu_l = liquid.dynamic_viscosity
    roughness = pipe.roughness
    grav = pipe.gravity
    beta = pipe.inclination

    # Set inital guess that solves from the top
    # the solution space is difficult
    # TODO find how to set a good initial guess
    height_initial = np.ones_like(u_gs) * 0.95
    # the critical heights at which waves would start to grow
    height_tilde = newton(
        equations.stratified.wave_growth, height_initial, args=(u_gs, liquid, gas, pipe)
    )
    # dimensional
    geom = Geometry(
        height_tilde, non_dimensional=False, pipe=pipe, u_gs=u_gs, u_ls=u_ls
    )

    # right hand side
    # get the actual fluid average reynolds and related frtiction factor
    reynolds_l_actual = rho_l * geom.vel_l * geom.hydr_diam_l / mu_l
    friction_l = friction_factor.niazkar_and_churchill(reynolds_l_actual, roughness)

    rhs = grav * pipe.diameter * (1 - height_tilde) * np.cos(beta) / friction_l

    # left hand side with actual fluid velocity
    lhs = geom.vel_l ** 2

    # the condition
    return lhs > rhs
