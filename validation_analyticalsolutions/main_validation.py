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

# iterate over the three different boundary conditions
for bc in boundary_conditions:
    
    print(f"Calculating for {bc} boundary condition")        
        
    if bc == "Dirichlet":
        
        # store temperature profiles for this boundary condition
        T_analytical[bc] = {}
        T_numerical[bc] = {}
        
        # constant surface temperature at T_surface = 800 K
        T_surface = 800
        
        # iterate over each time step to calculate and store the analytical and numerical temperatures
        for t in range(1,time_total):
                    
            # calculate the analytical solution for the temperature profile
            T_ansol = T_surface + (T_initial - T_surface)*special.erf(x_grid/(2 * np.sqrt(alpha * t)))
            # save temperature data into the analytical temperature
            T_analytical[bc][t] = T_ansol
            
            # calculating the numerical solution for the temperature profile
            T_num = 0
            # save temperature data into numerical solution
            T_numerical[bc][t] = T_num

     
    elif bc == "Neunman":
        
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

            # calculating the numerical solution for the temperature profile
            T_num = 0
            # save temperature data into numerical solution
            T_numerical[bc][t] = T_num
        
    elif bc == "Robin":
        
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
            
            # calculating the numerical solution for the temperature profile
            T_num = 100
            # save temperature data into numerical solution
            T_numerical[bc][t] = T_num

# condense all temperatures into one dictionary
validation_temperatures = {"analytical": T_analytical, "numerical": T_numerical}

# save in a pickle to retrieve and plot later
with open('validation_temperatures.pickle', 'wb') as handle:
    pickle.dump(validation_temperatures, handle)
    print("Temperature data saved into validation_temperatures.pickle")
    
grid_time = [x_grid, time_total]
with open('grid_time.pickle', 'wb') as handle:
    pickle.dump(grid_time, handle)