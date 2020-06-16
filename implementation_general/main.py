"""
This script uses the CN scheme to numerically solve the heat diffusion equation.

Boundary conditions
------------------
Surface: 
a) q_inc - q _conv - q_rad (where q_inc is a function of time)
b) q_inc - q_losses (using a total heat transfer coefficient)
Unexposed: 
a) conductive losses to an aluminium block
b) insulated boundary


1-D. Inet solid. Homogoneous. Constant properties.
"""

# import python libraries
import numpy as np
import pickle
import time

# import my own function
from cn import general_temperatures

# analyse three types of incident heat fluxes
heat_fluxes = ["Constant", "Linear", "Quadratic"]

# create list with boundary conditions to be validated
boundary_conditions_surface = ["Linear", "Non-linear"]
boundary_conditions_back = ["Insulation", "Aluminium-block"]

# create a dictionary to store temperature profiles for all boundary conditions and all times
Temperatures = {}

# parameters for the calculations
T_initial = 288                      # K
T_air = T_initial                    # K
time_total = 101                      # s
sample_length = 0.025                # m
space_divisions = 100                # -
h = 45                               # W/m2K for linearised surface bc with constant heat transfer coefficient

# define a range of properties to be used to evaluate the response of the material (these are syntethic properties)
k = np.linspace(0.2, 0.5, 4)           # W/mK
alpha = np.linspace(1e-7, 1.0e-6, 4)   # J/kgK

# create spatial mesh
dx = sample_length/(space_divisions - 1)
x_grid = np.array([i * dx for i in range(space_divisions)])

# define time step based on the spatial mesh for each alpha and create mesh
dt_all = (1 / 3) * (dx**2 / alpha)
time_divisions = time_total / dt_all
t_grid = []
for i,dt in enumerate(dt_all):
    t_grid.append(np.array([n * dt for n in range(int(time_divisions[i]))]))
    
upsilon = (alpha*dt)/(2*dx**2)

# iterate over three possible types of heat fluxes
for hf_type in heat_fluxes:
    
    print("---------")
    print(f" {hf_type} heat flux ")
    print("---------")
    Temperatures[f"hf-type:_{hf_type}"] = {}
    
    # according to the type of heat flux, different parameters are used to define the function
    q_all = []
    if hf_type == "Constant":
        q = np.linspace(20,80,4)
        q_all.append(q)
    elif hf_type == "Linear":
        q = np.linspace(0.1,0.7,4)
        q_all.append(q)
    elif hf_type == "Quadratic":
        q = np.linspace(2.5e-5, 1e-4,4)
        q_all.append(q)

    # iterate over surface boundary conditions
    for bc_surface in boundary_conditions_surface:
        
        print(f"-- Surface boundary condition {bc_surface} --")
        Temperatures[f"hf-type:_{hf_type}"][f"Surface_{bc_surface}"] = {}
        
        # iterate over back boundary conditions
        for bc_back in boundary_conditions_back:
                        
            start = time.time()
            Temperatures[f"hf-type:_{hf_type}"][f"Surface_{bc_surface}"][f"Back_{bc_back}"] = general_temperatures(hf_type, T_initial, T_air, time_total, k, alpha, dx, x_grid,
                        space_divisions,dt_all,t_grid, upsilon, bc_surface, bc_back, q, h)
            
            print(f"Time taken for {bc_surface} and {bc_back}: {np.round(time.time() - start,2)} seconds")


# condense all data to be saved including important plotting parameters
total_data = {"Temperatures": Temperatures, "extra_data": {
        "t_grid":t_grid, "x_grid":x_grid, "time_total": time_total, "alpha": alpha, "k": k, "q_all": q_all},
    }

# save in a pickle to retrieve and plot later
with open('total_data_general.pickle', 'wb') as handle:
    pickle.dump(total_data, handle)
    print("All data saved into total_data_general.pickle")
    