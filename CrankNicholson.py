# ------ 1D HEAT TRANSFER MODEL USING CRANK-NICHOLSON SCHEME. ACCOMODATES DIFFERENT BOUNDARY CONDITIONS -----------

# --- IMPORT PYTHON LIBRARIES

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb  # python debugger

# --- IMPORT MY OWN DEFINED FUNCTIONS

import myplotting as myplot
import CN_matrix as mymatrix
import mesh_icond as mymesh
import material_properties as mtp

# --- DEFINE SAMPLE CHARACTERISTICS

sample_depth = 0.025          # meters
space_mesh_divisions = 101    # (-)
time_duration = 900           # seconds
sample_material = "PMMA"      # PMMA, PA6, TIMBER, ALUMINIUM
initial_temperature = 293     # C
aluminium_depth = 0.02        # meters

# ---- DEFINE PLOTTING LIMITS

x_lim_inicond = [0, sample_depth]
y_lim_inicond = [0, initial_temperature + 25]
x_lim_other = [0, sample_depth]
y_lim_other = [0, 800]
figure_size = (16, 12)

# --- DEFINE THE BOUNDARY CONDITIONS

BC_surface = "conv_rad"                              # "const_temp", "const_nhf", "convection", "conv_rad"
BC_back = "insulated"                                 # "semi-inf", "insulated" or "al_block"
solving_method = "n_based"                           # "n_based" or "n+1_based". How surface radiation is considered for t = n+1
surface_temperature = 473                            # Define the surface temperature in C
surface_nhf = 20                                     # Define surface nhf in kW/m2
air_temperature = 293                                # K
h_convective = 0.015                                 # kW/m2K

# Dictionary contains arguments that vary depending on the chosen boundary condition
BC_values = {"const_temp": (surface_temperature, initial_temperature),
             "const_nhf": (surface_nhf, initial_temperature),
             "convection": (h_convective, initial_temperature, air_temperature),
             }

# ---- MAIN FUNCTION TO SOLVE FOR THE TEMPERATURE PROFILE


def main_solver():
    # 1. Create list to hold all temperature values
    Temperature = []

    # 2. Define mesh and initial conditions
    x_grid, t_grid, dx, dt, material, T, Tn, sigma, time_mesh_divisions = mymesh.define_mesh_icond(sample_depth, space_mesh_divisions, time_duration, sample_material, initial_temperature)
    Temperature.append(T)

    # 3. Define the Incident Heat Flux (IHF) (defined as an array so it requires definition of time_mesh_divisions)
    q_incident = material["absorptivity"] * 90 + np.zeros(t_grid.shape[0] + 1)  # kW/m2
    BC_values["conv_rad"] = (h_convective, initial_temperature, air_temperature, q_incident)

    # 4. Define matrix A
    A = mymatrix.tridiag_matrix(sigma, space_mesh_divisions, (BC_surface, BC_back), BC_values[BC_surface], material, dx, T, solving_method)
    print(A)
    # 5. Iterate to calculate T(x) for every time
    for i, _ in enumerate(t_grid):
        b = mymatrix.vector_b(sigma, space_mesh_divisions, T, (BC_surface, BC_back), BC_values[BC_surface], dx, material, i, solving_method)
        Tn = np.linalg.solve(A, b)
        Temperature.append(Tn)
        T = Tn.copy()
        print(T[0])
        print(i * dt)
        if np.isnan(b[0]):
            return None

    # Plot final temperature gradient
    myplot.plot_tempgrad(Tn, x_grid, figure_size, x_lim_other, y_lim_other, "Final Temperature Gradient", None, "show")

    # Plot verification plot if BC_surface = "const_temp", "const_nhf" or "convection" and B_back = "semi_inf"
    myplot.verify_plot(Tn, x_grid, x_lim_other, y_lim_other, time_duration, (BC_surface, BC_back), BC_values[BC_surface], material)

    return Temperature


# --- Call the solver and evaluate results
Temperature = main_solver()
