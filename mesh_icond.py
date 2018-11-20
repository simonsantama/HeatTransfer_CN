# ---- DEFINES THE MESH AND THE INITIAL CONDITIONS FOR THE PROBLEM

import material_properties as mtp
import numpy as np


def define_mesh_icond(sample_depth, space_mesh_divisions, time_duration, material, initial_temperature):
    """"
    Defines the grid, chooses material, initializes matrices.

    Returns space and time grids, material properties, T, Tn and sigma
    """
    dx = sample_depth / (space_mesh_divisions - 1)
    x_grid = np.array([i * dx for i in range(space_mesh_divisions)])

    # Material
    material = material.upper()
    if material == "PMMA":
        material = mtp.PMMA
    elif material == "PA6":
        material = mtp.PA6
    elif material == "TIMBER":
        material = mtp.TIMBER
    material["alpha"] = material["k"] / (material["rho"] * material["c"])

    # dt is defined following Von Neumann stability analysis
    dt = (1 / 3) * (dx**2 / material["alpha"])
    time_mesh_divisions = time_duration / dt
    t_grid = np.array([n * dt for n in range(int(time_mesh_divisions))])

    sigma = material["alpha"] * dt / (2 * dx**2)

    # Temperatures
    T = np.zeros(space_mesh_divisions) + initial_temperature
    Tn = np.empty_like(T)

    return (x_grid, t_grid, dx, dt, material, T, Tn, sigma, time_mesh_divisions)
