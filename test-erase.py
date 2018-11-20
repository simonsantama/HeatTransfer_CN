import numpy as np
import CN_matrix as mymatrix
import mesh_icond as mymesh
import pdb

sigma = 0.5
sample_depth = 0.025          # meters
space_mesh_divisions = 101    # (-)
time_duration = 600           # seconds
sample_material = "PMMA"      # PMMA, PA6, TIMBER, ALUMINIUM
initial_temperature = 293     # C
space_mesh_divisions = 101
BC_surface = "conv_rad"                              # "const_temp", "const_nhf", "convection", "conv_rad"
BC_back = "semi_inf"                                 # "semi-inf", "insulated" or "al_block"
surface_temperature = 473                            # Define the surface temperature in C
surface_nhf = 20                                     # Define surface nhf in kW/m2
air_temperature = 293                                # K
h_convective = 0.015                                 # kW/m2K

# Dictionary contains arguments that vary depending on the chosen boundary condition
BC_values = {"const_temp": (surface_temperature, initial_temperature),
             "const_nhf": (surface_nhf, initial_temperature),
             "convection": (h_convective, initial_temperature, air_temperature),
             }

x_grid, t_grid, dx, dt, material, T, Tn, sigma, time_mesh_divisions = mymesh.define_mesh_icond(sample_depth, space_mesh_divisions, time_duration, sample_material, initial_temperature)

# 3. Define the Incident Heat Flux (IHF) (defined as an array so it requires definition of time_mesh_divisions)
q_incident = 60 + np.zeros(t_grid.shape[0] + 1)  # kW/m2
BC_values["conv_rad"] = (h_convective, initial_temperature, air_temperature, q_incident)

for i in range(50):
    b = mymatrix.vector_b(sigma, space_mesh_divisions, T, (BC_surface, BC_back), BC_values[BC_surface], dx, material, i)
    T = T + 5000
print(T)
print(b)
