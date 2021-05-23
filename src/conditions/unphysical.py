"""
conditions which are impossible because of simple physics and not because
of any transition
"""
import ..general as general

def total_holdup(u_gs, u_ls, liquid, gas, pipe):
    """
    check locations which are unphysical because the total equivalent holdup 
    is larger than one
    """
    liquid_holdup = general.fluid_area_ratio(u_ls, liquid, pipe)
    gas_holdup = general.fluid_area_ratio(u_gs, gas, pipe)

    return liquid_holdup+gas_holdup>1