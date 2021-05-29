"""main module to use the other modules
"""

if __name__ == "__main__":
    """
    example of usage
    """
    import matplotlib.pyplot as plt
    import generate_data
    import fluids
    import parse_maps
    import visualization

    total_mass_flow = 2
    quality = 0.1
    liq_massflow = total_mass_flow * (1 - quality)
    gas_massflow = total_mass_flow * quality

    ugs_temp, uls_temp = generate_data.generate_velocity_maps()
    liq_temp = fluids.Liquid(
        density=998,
        bubble_surface_tension=0.073,
        mass_flowrate=liq_massflow,
        dynamic_viscosity=8.9e-4,
    )
    gas_temp = fluids.Gas(
        density=1.225, mass_flowrate=gas_massflow, dynamic_viscosity=18.3e-6
    )

    for inclination in [1]:  # [-90, -80, -30, -1, 0, 1, 30, 80, 90]:
        pipe_temp = fluids.Pipe(diameter=0.3, inclination=inclination, roughness=0.001)
        categories = parse_maps.get_categories_maps(
            ugs_temp, uls_temp, liq_temp, gas_temp, pipe_temp
        )

        fig_temp, ax = visualization.plot_map(
            categories,
            liq_temp,
            gas_temp,
            pipe_temp,
            x_ticks=ugs_temp[0, :],
            y_ticks=uls_temp[:, 0],
        )

    plt.show(block=True)
