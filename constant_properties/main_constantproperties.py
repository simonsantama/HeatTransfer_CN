"""
Main function for implementation of Crank-Nicolson scheme.
Constant thermal properties
"""

# import libraries
import numpy as np
import sys

# add paths to import the CN solver with insulation or with aluminium block as the back face boundary condition
sys.path.insert(1, r"C:\\Users\\s1475174\\Documents\\Python_Projects\\HeatTransfer_CrankNicholson\\HeatTransfer_CN\\constant_properties\\backface_aluminiumblock")
sys.path.insert(1, r"C:\\Users\\s1475174\\Documents\\Python_Projects\\HeatTransfer_CrankNicholson\\HeatTransfer_CN\\constant_properties\\backface_insulation")

from cn_insulation import general_solver as solver_insulation
from cn_aluminium import general_solver as solver_aluminium

def main_CN_constantproperties(heat_flux, surface_losses, backface_losses, time, space, properties, 
                               temperatures_initial):
    """
    This is the main function used to define the heat transfer scenario and call the corresponding algorithms to
    use the CN scheme to calculate the temperature evolution in the solid. 

    It accomodates different surface and back face boundary conditions, as well as multiple materials.

    Parameters
    ----------
    heat_flux : list
        Type of IHF and necessary parameters to calculate an array of IHF in time.
    surface_losses : list
        Type of surface losses to be considered and necessary parameters
    backface_losses: list
        Type of backface losses to be considered and necessary parameters
    time : int
        Total testing time
    space : np.array
        Spatial grid and necessary data for implemententing CN with two materials
    properties : list
        list of thermal properties needed
    temperatures_initial: list
        Initial solid phase temperature and gas phase temperture

    Returns
    -------
    all_data: dict
        temperature_profile, surface_temperature and nhf.

    """
    # create dictionary with all the data that is passed after all the calculations are completed
    all_data = {}
    
    # determine which solver to use depending on the backface boundary condition
    if backface_losses[0] == "insulated":
        solver = solver_insulation
    elif backface_losses[0] == "aluminium":
        solver = solver_aluminium
    
    # extract the spatial domain and create the spatial grid
    x_grid = np.linspace(0, space[0], space[1])
    all_data["x_grid"] = x_grid

    # extract initial and air temperatures
    temperature_initial = temperatures_initial[0]
    temperature_air = temperatures_initial[1]
    
    # calculate upsilon and thermal diffusivity alpha
    rho = properties[0]
    k = properties[1]
    c = properties[2]
    alpha = k/rho/c   
    all_data["properties"] = properties

    # calculate the temporal grid
    dx = x_grid[1] - x_grid[0]
    dt = (1 / 3) * (dx**2 / alpha)
    t_grid = np.arange(0,time[0],dt)
    all_data["t_grid"] = t_grid
    upsilon = (alpha*dt)/(2*dx**2)

    # filter the shape of the heat flux array according to the type of heat flux used
    hf = []
    if heat_flux[0] == "constant":
        for i in heat_flux[1]:
            hf.append(i + np.zeros_like(t_grid))
    elif heat_flux[0] == "linear":
        for i in heat_flux[1]:
            hf.append(i*t_grid)
    elif heat_flux[0] == "quadratic":
        for i in heat_flux[1]:
            hf.append(i*t_grid**2)
    elif heat_flux[0] == "cubic":
        pass
    
    all_data["heat_fluxes"] = hf
    
    # iterate over every heat flux given
    for i, hf_array in enumerate(hf):
        
        # calculate temperature profiles in time
        calculated_temperatures = solver(heat_flux = hf_array*1000, temp_initial = temperature_initial, 
                                          temp_air = temperature_air, k = k, alpha = alpha,
                                          x_grid = x_grid, t_grid = t_grid, upsilon = upsilon,
                                          bc_surface = surface_losses, sigma = 5.67e-8 )
        
        all_data[f"IHF_{heat_flux[0]}_{heat_flux[1][i]}"] = calculated_temperatures
    
    return all_data