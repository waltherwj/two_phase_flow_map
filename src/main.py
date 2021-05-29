"""main module to use the other modules
"""

if __name__ == "__main__":
    import numpy as np
    from generate_data import generate_velocity_maps
    from fluids import Liquid, Gas, Pipe
    from parse_maps import get_categories_maps

    total_mass_flow = 0.01
    quality = 0.1
    liq_massflow = total_mass_flow * (1 - quality)
    gas_massflow = total_mass_flow * quality

    ugs_temp, uls_temp = generate_velocity_maps()
    liq_temp = Liquid(
        density=998,
        bubble_surface_tension=0.073,
        mass_flowrate=liq_massflow,
        dynamic_viscosity=8.9e-4,
    )
    gas_temp = Gas(density=1.225, mass_flowrate=gas_massflow, dynamic_viscosity=18.3e-6)

    for inclination in [-90, -80, -30, -1, 0, 1, 30, 80, 90]:
        pipe_temp = Pipe(diameter=5.1e-2, inclination=inclination, roughness=0.00001)
        categories = get_categories_maps(
            ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
        )

        fig_temp, ax = plot_map(
            categories,
            liq_temp,
            gas_temp,
            pipe_temp,
            x_ticks=ugs_temp[0, :],
            y_ticks=uls_temp[:, 0],
        )

    plt.show(block=True)
