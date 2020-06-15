"""
This script uses the CN scheme to numerically solve the heat diffusion equation and compares the results to analytical
solutions for three different boundary conditions.

Analytical solutions are taken from Fundamentals of Heat and Mass Transfer (Incropera, Dewitt)

It assumes a 1 dimensional, semi-infinite, inert and homogenous solid with constant properties.
"""

# import python libraries
import numpy as np
from scipy import special
import pickle
import time

# import my own function
from cn_validation import cn_solver

# create list with boundary conditions to be validated
boundary_conditions = ["Dirichlet", "Neunman", "Robin"]

# create a dictionary to store temperature profiles for all boundary conditions and all times
T_analytical = {}
T_numerical = {}

# parameters for the calculations
T_initial = 300          # K
T_air = T_initial        # K
time_total = 901         # s
sample_length = 0.2      # m
space_divisions = 1000   # -

alpha = 1e-6            # m2/s - random value (close to PMMA)
k = 2e-1                # W/mK - random value (close to PMMA)

# create mesh (space) for the numerical solution
dx = sample_length / (space_divisions - 1)
x_grid = np.array([i * dx for i in range(space_divisions)])

# define dt according to Von Neunman stability analysis
dt = (1 / 3) * (dx**2 / alpha)
time_divisions = time_total / dt
t_grid = np.array([n * dt for n in range(int(time_divisions))])

upsilon = (alpha*dt)/(2*dx**2)

# iterate over the three different boundary conditions
for bc in boundary_conditions:
    
    print(f"Calculating for {bc} boundary condition")        
        
    # constant surface temperature
    if bc == "Dirichlet":
        
        start_dirich_an = time.time()
        # store temperature profiles for this boundary condition
        T_analytical[bc] = {}
        
        # constant surface temperature at T_surface = 800 K
        T_surface = 800
        
        # calculate analytical solution
        for t in range(1,time_total):
                    
            # calculate the analytical solution for the temperature profile
            T_ansol = T_surface + (T_initial - T_surface)*special.erf(x_grid/(2 * np.sqrt(alpha * t)))
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol
        print(f"Time taken for {bc} analytical solution: {np.round(time.time() - start_dirich_an,2)} seconds")
        
        
        start_dirich_num = time.time()
        # calculate the numerical solution
        bc_data = T_surface
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data, space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num
        print(f"Time taken for {bc} numerical solution: {np.round(time.time() - start_dirich_num,2)} seconds")
     
    # constant surface heat flux
    elif bc == "Neunman":
        
        start_neunm_an = time.time()
        # store temperature profiles for this boundary condition
        T_analytical[bc] = {}
        T_numerical[bc] = {}
        
        # constant surface heat flux of q = 5 kW/m2
        q = 5000
        
        # iterate over each time step to calculate and store the analytical and numerical temperatures
        for t in range(1,time_total):
                    
            # calculate the analytical solution for the temperature profile
            T_ansol = T_initial + \
                    (2*q/k)*np.sqrt(alpha*t/np.pi)*np.exp(-x_grid**2/(4 * alpha * t)) - \
                    (q*x_grid/k)*special.erfc(x_grid/(2 * np.sqrt(alpha * t)))
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol
        print(f"Time taken for {bc} analytical solution: {np.round(time.time() - start_neunm_an,2)} seconds")

        start_neum_num = time.time()
        # calculate the numerical solution
        bc_data = q
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data, space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num
        print(f"Time taken for {bc} numerical solution: {np.round(time.time() - start_neum_num,2)} seconds")
        
    # surface convective heating
    elif bc == "Robin":
        
        start_robin_an = time.time()
        # store temperature profiles for this boundary condition
        T_analytical[bc] = {}
        T_numerical[bc] = {}
        
        # surface convection with a convective heat transfer coefficient of 50 W/m2K and a gas temperature of 800 K
        h = 50               # W/m2K
        T_infinity = 800     # K
        
        # iterate over each time step to calculate and store the analytical and numerical temperatures
        for t in range(1,time_total):
                    
            # calculate the analytical solution for the temperature profile
            T_ansol = T_initial + \
                    (T_infinity - T_initial)*(\
                    special.erfc(x_grid/(2*np.sqrt(alpha*t))) - \
                    np.exp(h*x_grid/k + h**2*alpha*t/k**2)*special.erfc(x_grid/(2*np.sqrt(alpha*t)) + h*np.sqrt(alpha*t)/k) \
                    )

            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol    
        print(f"Time taken for {bc} analytical solution: {np.round(time.time() - start_robin_an,2)} seconds")
        
        start_robin_num = time.time()
        # calculate the numerical solution
        bc_data = h, T_infinity
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data, space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num
        print(f"Time taken for {bc} numerical solution: {np.round(time.time() - start_robin_num,2)} seconds")

# since analytical solutions are evaluted at 1 Hz, drop columns from numerical solution to facilitate plotting
T_numerical_int = {}
list_indexes = []

for bc in boundary_conditions:
    T_numerical_int[bc] = {}
    for i in range(1,time_total):
        for key in T_numerical[bc].keys():
            if float(key) > i:
                T_numerical_int[bc][i] = T_numerical[bc][key]
                break
            continue 
    
# condense all temperatures into one dictionary
validation_temperatures = {"analytical": T_analytical, "numerical": T_numerical_int}

# save in a pickle to retrieve and plot later
with open('validation_temperatures.pickle', 'wb') as handle:
    pickle.dump(validation_temperatures, handle)
    print("Temperature data saved into validation_temperatures.pickle")
    
grid_time = [x_grid, time_total]
with open('grid_time.pickle', 'wb') as handle:
    pickle.dump(grid_time, handle)