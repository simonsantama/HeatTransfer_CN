import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

# 1-D Heat Transfer Model. Uses Crank-Nicholson scheme

# --- 1. Define material properties ---
PMMA = {"k": 0.0002,  # kW/mK [Vermesi_PhD]
        "rho": 1190,  # kg/m3 [Vermesi_PhD]
        "c": 1.606,  # kJ/kgK [Vermesi_PhD]
        "emissivity": 1,
        "absortivity": 1,
        }

PA6 = {"k": 1,
       "rho": 1,
       "c": 1,
       "emissivity": 1,
       "absortivity": 1,
       }

TIMBER = {"k": 1,
          "rho": 1,
          "c": 1,
          "emissivity": 1,
          "absortivity": 1,
          }

ALUMINIUM = {"k": 0.167,  # kW/mK [Vermesi_PhD]
             "rho": 2700,  # kg/m3 [Vermesi_PhD]
             "c": 0.896,  # kJ/kgK [Vermesi_PhD]
             "m": 0.4,    # kg
             "emissivity": 1,
             "absortivity": 1,
             }

# --- 2. Define sample characteristics ---
sample_depth0 = 0.025  # meters
space_mesh_points0 = 101    # (-)
duration_experiment0 = 900    # seconds
time_mesh_points0 = 1000   # (-)
sample_material0 = "PMMA"  # PMMA, PA6, TIMBER, ALUMINIUM
initial_temperature0 = 20     # C

# Aluminium characteristics (assume lumped capacitance)
sample_depth1 = 0.02                   # meters
sample_material1 = "ALUMINIUM"            # ALUMINIUM
initial_temperature1 = initial_temperature0   # C

# --- 3. Define the boundary conditions ---
BC_type = "const_temp"    # "const_temp", "const_nhf", "convection", "conv_rad"
BC_type_back = "semi_inf"      # "semi-inf" or "al_block"

if BC_type == "const_temp" and BC_type_back == "semi_inf":
    t_surface = 200  # C
    conditions_for_b = (t_surface, initial_temperature0)
elif BC_type == "const_temp":
    nhf_surface = 20  # kW/m2
else:
    q_inc = 25 + np.zeros(time_mesh_points0)  # np.array(time_mesh_points)
    temperature_gas = 20 + np.zeros(time_mesh_points0)  # np.array(time_mesh_points)

# --- 4. Some plotting characteristics ---
mpl.rc('font', size=18)
mpl.rc('font', family='Arial')
figure_size = (16, 12)
x_lim_inicond = [0, sample_depth0]
y_lim_inicond = [0, initial_temperature0 + 25]
x_lim_other = [0, sample_depth0]
y_lim_other = [0, 550]

# --- 5. Define the mesh and initial condition ---


def define_mesh_icond(length_meters, length_divisions, time_domain, time_divisions, material, initial_temperature):
    """"
    Defines the grid, chooses material, initializes matrices.

    Returns space and time grids, material properties, T, Tn and sigma
    """
    dx = length_meters / (length_divisions - 1)
    dt = time_domain / (time_divisions - 1)
    x_grid = np.array([i * dx for i in range(length_divisions)])
    t_grid = np.array([n * dt for n in range(time_divisions)])
    # Material
    material = material.upper()
    if material == "PMMA":
        material = PMMA
    elif material == "PA6":
        material = PA6
    elif material == "TIMBER":
        material = TIMBER
    material["alpha"] = material["k"] / (material["rho"] * material["c"])
    sigma = material["alpha"] * dt / (2 * dx**2)
    # Temperatures
    T = np.zeros(length_divisions) + initial_temperature
    Tn = np.empty_like(T)

    return (x_grid, t_grid, dx, dt, material, T, Tn, sigma)

# --- 6. Define function to plot the temperature gradient ---


def plot_tempgrad(T, x_grid, figure_size, x_lim, y_lim, title, save_format):
    """
    Plots temperature gradient at a given time

    Save format = None for no saving
    """
    fig, ax = plt.subplots(figsize=figure_size)
    ax.plot(x_grid, T, linewidth=2.5)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_xlabel("Depth [m]")
    ax.set_ylabel("Temperature [C]")
    ax.set_title(title)
    if save_format == None:
        pass
    else:
        fig.savefig(title + save_format, dpi=300)

    return None
