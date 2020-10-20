"""
Solver for the Crank-Nicolson implementation of heat transfer
One dimension
Constant properties
Insulated back face
Accomodates linear or non-linear surface boundary condition
"""

# import libraries
import numpy as np
import pandas as pd

def general_solver(heat_flux, temp_initial, temp_air, k, alpha, x_grid,t_grid, upsilon, 
                   bc_surface, sigma):
    """
    General function to implement CN.
    Returns dictionary that includes all calculated temperatures for given boundary condition
    
    CN Scheme:
    ---------
        [A]{T^n+1} = [B]{T^n} = [b]
    
    
    Parameters:
    ----------
    heat_flux: np.array
        array of IHF in time
    
    temp_initial: int
        initial temperature in K
        
    temp_air: int
        air (infinity) temperature for calculating convective losses       
        
    k: int
        thermal conductivity in W/mK
        
    alpha: float
        thermal diffusivity in m2/s
    
    x_grid: np.array
        one dimentional spatial array

    t_grid: np.array
        array over which the algorithm advances in time
        
        
    upsilon: Fourier number divided by 2
        list
        
    bc_surface: list
        Type of surface losses to be considered and necessary parameters
        
    sigma: Stefan Boltzman constant
        float

    Returns:
    -------
    Temperature: pd.DataFrame
        one-dimensional temperature profile in time
    
    """

    # temperatures are reported as a data frame, where each column is a step in time
    temperatures = pd.DataFrame(columns = [n for n in t_grid])

    # extract the necessary parameters to determine the surface heat losses
    if bc_surface[0] == "linear":
        h = bc_surface[1] + bc_surface[2]
        hc = 0
        emissivity = 0
    elif bc_surface[0] == "non-linear":
        h = 0
        hc = bc_surface[1]
        emissivity = bc_surface[2]

    # initialize temperature arrays for present and future temperatures
    T = np.zeros_like(x_grid) + temp_initial
    Tn = np.zeros_like(x_grid)

    # iterate over each time step
    temperatures.iloc[:,0] = T
    for j, t in enumerate(t_grid[:-1]):
        
        # create tri-diagonal matrix A
        A = tridiag_matrix(bc_surface_type = bc_surface[0], upsilon = upsilon, 
                           space_divisions = len(x_grid), dx = x_grid[1] - x_grid[0], 
                           k = k, T = T, h = h, hc = hc, emissivity = emissivity, sigma = sigma)
        
        # create vector b
        b = vector_b(bc_surface_type = bc_surface[0], upsilon = upsilon, 
                     space_divisions = len(x_grid), dx = x_grid[1] - x_grid[0], 
                     k = k, T = T, T_air = temp_air, heat_flux = heat_flux, h = h, hc = hc, 
                              emissivity = emissivity, sigma = sigma, j = j)
        
        # calculate value of future temperature
        Tn = np.linalg.solve(A,b)
        
        # update present temperature
        T = Tn.copy()
        
        # store temperature profile at this time in the data frame
        temperatures.iloc[:, j+1] = Tn
            
    return temperatures   
    
    

# function to create tri-diagonal matrix
def tridiag_matrix(bc_surface_type, upsilon, space_divisions, dx, k, T, h, hc, emissivity, sigma):
    """
    Creates tridiagonal matrix A
    Linear system to be solved is Ax = b, and x represents temperature values at time n+1

    Parameters:
    ----------
    bc_surface_type: str
        boundary condition at the surface
    
    upsilon: float
        Fourier number divided by 2. Upsilon = alpha*dt/2*dx2
        
    space_divisions: int
        number of nodes in the spatial domain
        
    dx: float
        size of cell in space domain in m
        
    k: float
        thermal conductivity in W/mK
        
    T: np.array
        array of present temperatures
        
    h: int
        total heat transfer coefficient for the linearised surface boundary condition
        
    hc: int
        convective heat transfer coefficient
        
    emissivity: float
        surface emmisivity, assumed constant
        
    sigma: float
        Stefan Boltzman constant
    
    Return:
    ------
    
    A: matrix to be inverted
        np.array
    
    """
    # create tri-diagonal matrix
    A = np.diagflat([-upsilon for i in range(space_divisions - 1)], -1) +\
        np.diagflat([1 + 2 * upsilon for i in range(space_divisions)]) +\
        np.diagflat([-upsilon for i in range(space_divisions - 1)], 1)

    # adjust matrix depending on the boundary condition at the exposed surface
    if bc_surface_type == "linear":
        A[0,0] = 1 + 2*upsilon + 2*upsilon*dx*h/k
        A[0,1] = -2*upsilon
    
    elif bc_surface_type == "non-linear":
        A[0,0] = 1 + 2*upsilon + 2*dx*hc*upsilon/k+ 8*emissivity*sigma*dx*upsilon*T[0]**3/k
        A[0,1] = -2*upsilon
    
    # adjust matrix for the back boundary conditions
    A[-1, -2] = - 2 * upsilon
    A[-1, -1] = 1 + 2 * upsilon

    return A



def vector_b(bc_surface_type, upsilon, space_divisions, dx, k, T, T_air, heat_flux, h, hc, emissivity, sigma, j):
    """
    Calculates vector b. Right hand side of linear system of equations

    Parameters:
    ----------
    bc_surface_type: str
        boundary condition at the surface
    
    upsilon: float
        Fourier number divided by 2. Upsilon = alpha*dt/2*dx2
        
    space_divisions: int
        number of nodes in the spatial domain
        
    dx: float
        size of cell in space domain in m
        
    k: float
        thermal conductivity in W/mK
        
    T: np.array
        array of present temperatures
        
    T_air: int
        air (infinity) temperature for calculating convective losses 
        
    heat_flux: np.array
        full array of IHF vs time
        
    h: int
        total heat transfer coefficient for the linearised surface boundary condition
        
    hc: int
        convective heat transfer coefficient

    emissivity: float
        surface emmisivity, assumed constant
        
    sigma: float
        Stefan Boltzman constant
        
    j: int
        present iteration number
    
    Returns:
    -------
    b: np.array
        vector to solve linear system of equations
    """
    
    # matrix B, similar to matrix A but multiplies T at present
    B = np.diagflat([upsilon for i in range(space_divisions - 1)], -1) +\
        np.diagflat([1 - 2 * upsilon for i in range(space_divisions)]) +\
        np.diagflat([upsilon for i in range(space_divisions - 1)], 1)

    # Calculate vector b (matrix B dot product with present temperatures)
    b = np.zeros(space_divisions)
    b[1:-1] = B[1:-1, :].dot(T)
    
    # adjust vector for the front boundary condition
    if bc_surface_type == "linear":
        b[0] = 2*upsilon*T[1] + (1 - 2*upsilon - upsilon*2*dx*h/k)*T[0] + 4*upsilon*dx*h*T_air/k + \
            2*dx*upsilon/k * (heat_flux[j+1]+heat_flux[j])
    
    elif bc_surface_type == "non-linear":
        b[0] = 2*upsilon*T[1] + (1- 2*upsilon - 2*dx*hc*upsilon/k)*T[0] + 4*dx*hc*upsilon*T_air/k + \
            4*emissivity*sigma*dx*upsilon*T[0]**4/k + 2*dx*upsilon/k * (heat_flux[j+1]+heat_flux[j])
    
    # adjust vector for the back boundary condition
    b[-1] = (1 - 2*upsilon)*T[-1] + 2*upsilon*T[-2]

    return b    
    



