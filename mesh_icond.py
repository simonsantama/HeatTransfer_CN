# ---- DEFINES THE MESH AND THE INITIAL CONDITIONS FOR THE PROBLEM

import material_properties as mtp
import numpy as np


def define_mesh_icond(sample_depth, space_mesh_division, time_duration, time_mesh_divisions, material, initial_temperature):
    """"
    Defines the grid, chooses material, initializes matrices.

    Returns space and time grids, material properties, T, Tn and sigma
    """
    dx = sample_depth / (space_mesh_division - 1)
    dt = time_duration / (time_mesh_divisions - 1)
    x_grid = np.array([i * dx for i in range(space_mesh_division)])
    t_grid = np.array([n * dt for n in range(time_mesh_divisions)])

    # Material
    material = material.upper()
    if material == "PMMA":
        material = mtp.PMMA
    elif material == "PA6":
        material = mtp.PA6
    elif material == "TIMBER":
        material = mtp.TIMBER
    material["alpha"] = material["k"] / (material["rho"] * material["c"])
    sigma = material["alpha"] * dt / (2 * dx**2)

    # Temperatures
    T = np.zeros(space_mesh_division) + initial_temperature
    Tn = np.empty_like(T)

    return (x_grid, t_grid, dx, dt, material, T, Tn, sigma)