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
time_total = 601         # s
sample_length = 0.2      # m
space_divisions = 400    # -
alpha = 1e-6             # m2/s - random value (close to PMMA)
k = 2e-1                 # W/mK - random value (close to PMMA)

# values used for analytical solutions
T_surface = 800          # K constant surface temperature for Dirichlet boundary condition
q = 5000                 # W/m2 constant surface heat flux for Neunman boundary condition
h = 50                   # W/m2K convective heat transfer coefficient for Robin boundary condition
T_infinity = 800         # K gas temperature for convective heat losses with Robin boundary condition

bc_data = [T_surface, q, (h,T_infinity)]

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
        start = time.time()
        
        # -- analytical solution
        T_analytical[bc] = {}
        # calculate analytical temperature and save into dictionary
        for t in range(1,time_total):
            T_ansol = T_surface + (T_initial - T_surface)*special.erf(x_grid/(2 * np.sqrt(alpha * t)))
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol            
        
        # -- numerical solution
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data[0], space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num     
        
        print(f" - time taken for {bc} boundary condition: {np.round(time.time() - start,2)} seconds")
        

    # constant surface heat flux
    elif bc == "Neunman":
        start = time.time()

        # -- analytical solution
        T_analytical[bc] = {}
        # calculate analytical temperature and save into dictionary
        for t in range(1,time_total):                    
            T_ansol = T_initial + \
                    (2*q/k)*np.sqrt(alpha*t/np.pi)*np.exp(-x_grid**2/(4 * alpha * t)) - \
                    (q*x_grid/k)*special.erfc(x_grid/(2 * np.sqrt(alpha * t)))
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol

        # -- numerical solution
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data[1], space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num
        
        print(f" - time taken for {bc} boundary condition: {np.round(time.time() - start,2)} seconds")
        
    # surface convective heating
    elif bc == "Robin":
        
        start = time.time()
        # -- analytical solution
        T_analytical[bc] = {}
        # calculate analytical temperature and save into dictionary
        for t in range(1,time_total):
            T_ansol = T_initial + \
                    (T_infinity - T_initial)*(\
                    special.erfc(x_grid/(2*np.sqrt(alpha*t))) - \
                    np.exp(h*x_grid/k + h**2*alpha*t/k**2)*special.erfc(x_grid/(2*np.sqrt(alpha*t)) + h*np.sqrt(alpha*t)/k) \
                    )
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol  

        # -- numerical solution
        T_num = cn_solver(x_grid, t_grid, upsilon, bc, bc_data[2], space_divisions, T_initial, dx, k)
        # save temperature data into numerical solution
        T_numerical[bc] = T_num

        print(f" - time taken for {bc} boundary condition: {np.round(time.time() - start,2)} seconds")

# obtain numerical data at 1 Hz to facilitate plotting
T_numerical_1Hz = {}
list_indexes = []

for bc in boundary_conditions:
    T_numerical_1Hz[bc] = {}
    time_counter = 0
    
    times = list(T_numerical[bc].keys())
    times.sort()
    
    for t in times:
        if t > time_counter:
            T_numerical_1Hz[bc][time_counter] = T_numerical[bc][t]
            time_counter += 1
  
# condense all temperatures into one dictionary
validation_temperatures = {"analytical": T_analytical, "numerical": T_numerical_1Hz}

# save in a pickle to retrieve and plot later
with open('validation_temperatures.pickle', 'wb') as handle:
    pickle.dump(validation_temperatures, handle)
    print("Temperature data saved into validation_temperatures.pickle")
    
grid_time = [x_grid, time_total]
with open('grid_time.pickle', 'wb') as handle:
    pickle.dump(grid_time, handle)