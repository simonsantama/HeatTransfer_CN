import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys

# 1-D Heat Transfer Model. Uses Crank-Nicholson scheme

# --- 1. Define material properties ---


PMMA = {"k": 0.0002,  # kW/mK [Vermesi_PhD]
        "rho": 1190,  # kg/m3 [Vermesi_PhD]
        "c": 1.606,   # kJ/kgK [Vermesi_PhD]
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


sample_depth = 0.025          # meters
space_mesh_divisions = 101    # (-)
time_duration = 900           # seconds
time_mesh_divisions = 1000    # (-)
sample_material = "PMMA"      # PMMA, PA6, TIMBER, ALUMINIUM
initial_temperature = 20      # C

# Aluminium characteristics (assume lumped capacitance)
aluminium_depth = 0.02        # meters

# --- 3. Define the boundary conditions ---


BC_surface = "const_temp"         # "const_temp", "const_nhf", "convection", "conv_rad"
BC_back = "semi_inf"      # "semi-inf" or "al_block"

# Surface BC 1: Constant surface temperature
if BC_surface == "const_temp":
    t_surface = 200  # C
    if BC_back == "semi_inf":
        BC_tuple = (t_surface, initial_temperature)
    elif BC_back == "insulated":
        BC_tuple = (t_surface, 0)
    elif BC_back == "al_block":
        pass
# Surface BC 2: Constant NHF
elif BC_surface == "const_nhf":
    pass

# --- 4. Some plotting characteristics ---


mpl.rc('font', size=18)
mpl.rc('font', family='Arial')
figure_size = (16, 12)
x_lim_inicond = [0, sample_depth]
y_lim_inicond = [0, initial_temperature + 25]
x_lim_other = [0, sample_depth]
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


def plot_tempgrad(T, x_grid, figure_size, x_lim, y_lim, title, save_format, show):
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
    # Decide whether or not to show the graph
    if show == "show":
        plt.show()
    elif show == "not-show":
        pass
    # Decide whether or not to save the graph
    if save_format == None:
        pass
    else:
        fig.savefig(title + save_format, dpi=300)

    return None

# --- 7. Define Matrix A ---


def tridiag_matrix(sigma, length_divisions, BC_type):
    """
    Creates tridiagonal matrix A

    BC_type = (BC_type, BC_back)

    BC_type: "const_temp", "const_nhf", "convection", "conv_rad"
    BC_back: "semi_inf" or "al_block"
    """

    A = np.diagflat([-sigma for i in range(length_divisions - 1)], -1) +\
        np.diagflat([1 + 2 * sigma for i in range(length_divisions)]) +\
        np.diagflat([-sigma for i in range(length_divisions - 1)], 1)

    # Boundary conditions
    if BC_type[0] == "const_temp":
        A[0, 0] = 1
        A[0, 1] = 0
        if BC_type[1] == "semi_inf":
            A[-1, -2] = 0
            A[-1, -1] = 1
    return A

# --- 8. Define the vector b (Matrix B x T)


def vector_b(sigma, length_divisions, T, BC_type, BC_tuple):
    """
    Calculates vector at the right hand side of the algebraic equation

    BC_type = (BC_type, BC_back)

    BC_type: "const_temp", "const_nhf", "convection", "conv_rad"
    BC_back: "semi_inf" or "al_block"
    """
    B = np.diagflat([sigma for i in range(length_divisions - 1)], -1) +\
        np.diagflat([1 - 2 * sigma for i in range(length_divisions)]) +\
        np.diagflat([sigma for i in range(length_divisions - 1)], 1)

   # Calculate vector b
    b = B.dot(T)

    # Boundary conditions
    if BC_type[0] == "const_temp":
        b[0] = BC_tuple[0]
        if BC_type[1] == "semi_inf":
            b[-1] = BC_tuple[1]
    return b

# --- . Define the main function which calls all other functions define here


def main_solver():
    # Create Temperature list. Each entry corresponds to temperature profile for a given time.
    Temperature = []

    # Define mesh and initial conditions
    x_grid, t_grid, dx, dt, material, T, Tn, sigma = define_mesh_icond(sample_depth, space_mesh_divisions, time_duration, time_mesh_divisions, sample_material, initial_temperature)
    Temperature.append(T)

    # Plot the initial condition
    plot_tempgrad(T, x_grid, figure_size, x_lim_inicond, y_lim_inicond, "Initial condition", None, "not-show")

    # Define matrix A (CN discretization + Boundary conditions)
    A = tridiag_matrix(sigma, space_mesh_divisions, (BC_surface, BC_back))

    # Define the vector b (B x T)
    for _ in range(1000):
        b = vector_b(sigma, space_mesh_divisions, T, (BC_surface, BC_back), BC_tuple)
        Tn = np.linalg.solve(A, b)
        Temperature.append(Tn)
        T = Tn.copy()
    return Temperature


# --- Call the solver and evaluate results
Temperature = main_solver()
print(Temperature)
