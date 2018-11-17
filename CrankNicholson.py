# ------ 1D HEAT TRANSFER MODEL USING CRANK-NICHOLSON SCHEME. ACCOMODATES DIFFERENT BOUNDARY CONDITIONS -----------

# --- IMPORT PYTHON LIBRARIES

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

# --- IMPORT MY OWN DEFINED FUNCTIONS

import myplotting as myplot
import CN_matrix as mymatrix
import mesh_icond as mymesh
import material_properties as mtp

# --- DEFINE SAMPLE CHARACTERISTICS

sample_depth = 0.025          # meters
space_mesh_divisions = 101    # (-)
time_duration = 600           # seconds
time_mesh_divisions = 1000    # (-)
sample_material = "PMMA"      # PMMA, PA6, TIMBER, ALUMINIUM
initial_temperature = 20      # C
aluminium_depth = 0.02        # meters

# ---- DEFINE PLOTTING LIMITS

x_lim_inicond = [0, sample_depth]
y_lim_inicond = [0, initial_temperature + 25]
x_lim_other = [0, sample_depth]
y_lim_other = [0, 550]


# --- DEFINE THE BOUNDARY CONDITIONS

BC_surface = "convection"   # "const_temp", "const_nhf", "convection", "conv_rad"
BC_back = "semi_inf"        # "semi-inf" or "al_block"
surface_temperature = 200   # Define the surface temperature in C
surface_nhf = 20            # Define surface nhf in kW/m2
air_temperature = 450       # C
h_convective = 0.015        # kW/m2K
BC_values = {"const_temp": (surface_temperature, initial_temperature),
             "const_nhf": (surface_nhf, initial_temperature),
             "convection": (h_convective, initial_temperature, air_temperature)
             }


def main_solver():
    # 1. Create list to hold all temperature values
    Temperature = []

    # 2. Define mesh and initial conditions
    x_grid, t_grid, dx, dt, material, T, Tn, sigma = mymesh.define_mesh_icond(sample_depth, space_mesh_divisions, time_duration, time_mesh_divisions, sample_material, initial_temperature)
    Temperature.append(T)

    # 3. Plot the initial condition
    # myplot.plot_tempgrad(T, x_grid, figure_size, x_lim_inicond, y_lim_inicond, "Initial condition", None, "not-show")

    # 4. Define matrix A
    A = mymatrix.tridiag_matrix(sigma, space_mesh_divisions, (BC_surface, BC_back), BC_values[BC_surface], material, dx)

    # 5. Iterate to calculate T(x) for every time
    for _ in t_grid:
        b = mymatrix.vector_b(sigma, space_mesh_divisions, T, (BC_surface, BC_back), BC_values[BC_surface], dx, material)
        Tn = np.linalg.solve(A, b)
        Temperature.append(Tn)
        T = Tn.copy()

    # Plot final temperature gradient
    # plot_tempgrad(Tn, x_grid, figure_size, x_lim_other, y_lim_other, "Final Temperature Gradient", "pdf", "not-show")

    # Plot verification plot if BC_surface = "const_temp", "const_nhf" or "convection" and B_back = "semi_inf"
    myplot.verify_plot(Tn, x_grid, x_lim_other, y_lim_other, time_duration, (BC_surface, BC_back), BC_values[BC_surface], material)

    return Temperature


# --- Call the solver and evaluate results
# Temperature = main_solver()

for time_duration in [200, 600, 800]:
    main_solver()
    print(time_duration)
