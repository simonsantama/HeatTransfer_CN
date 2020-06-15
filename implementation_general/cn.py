"""
Crank-Nicholson implementation of 1D heat transfer
"""

# import libraries
import numpy as np


def general_temperatures(hf_type, T_initial, T_air, time_total, k, alpha, dx, x_grid, space_divisions, dt_all,
                    t_grid, upsilon, bc_surface, bc_back):
    """
    General function to implement CN.
    Creates different functions of time for the incident heat flux.
    Returns dictionary that includes all calculated temperature for given boundary conditions evaluated at different
    properties and parameters
    
    Validation function: Exposed surface uses dirichlet, neunman or robin boundary conditions.
    Semi-infinite solid
    
    CN Scheme:
    ---------
        [A]{T^n+1} = [B]{T^n} = [b]
    
    
    Parameters:
    ----------
    hf_type: type of heat flux function for the incident heat flux.
        str
    
    T_initial: initial temperature in K
        int
        
    T_air: air(infinity) temperature in K for convective losses. Usually equals T_initial but not necessary
        int        

    time_total: total time for calculations
        int
        
    k: thermal conductivity in W/mK. Array with range of values
        np.array
        
    alpha: thermal diffusivity in m2/s. Array with range of values
        np.array
        
    dx: size of cell in space domain in m
        float 
    
    x_grid: one dimensional spatial domain
        np.array

    space_divisions: number of nodes in the spatial domain
        int

    dt_all: different thermal diffusivities define different time steps for a given dx
        np.array

    t_grid: list where each entry is the temporal domain for a given alpha and dt
        list
        
    upsilon: Fourier number divided by 2
        float
        
    bc_surface: surface boundary condition (linear or non-linear)
        str
        
    bc_back: back boundary condition (insulation or conductive losses to aluminium block)
        str

    Returns:
    -------
    Temperature: dictionary with temperature arrays for different alpha and k
        dict
    
    """

    # according to the type of heat flux, different parameters are used to define the function
    if hf_type == "Constant":
        q = np.linspace(20,80,4)
    elif hf_type == "Linear":
        q = np.linspace(100,700,4)
    elif hf_type == "Quadratic":
        q = np.linspace(2.5e-5, 1e-4,4)

    temperatures = {}

    # iterate over the different heat flux defined above
    for heat_flux in q:
    
        temperatures[f"q:_{heat_flux}"] = {}
        
        # iterate over all value pairs of alpha and k
        for i in range(len(alpha)):
            k_this = k[i]
            upsilon_this = upsilon[i]
            alpha_this = alpha[i]
            t_grid_this = t_grid[i]
            
            # create tridiagonal matrix A
            A = tridiag_matrix(bc_surface, bc_back, upsilon_this, space_divisions, dx, k_this)
            
            # initialise temperature arrays for present and future temperatures
            T = np.zeros_like(x_grid) + T_initial
            Tn = np.zeros_like(x_grid)
            
            # define the incident heat flux
            if hf_type == "Constant":
                q_array = np.zeros_like(t_grid_this) + heat_flux
            elif hf_type == "Linear":
                q_array = heat_flux * t_grid_this
            elif hf_type == "Quadratic":
                q_array = heat_flux * t_grid_this * t_grid_this
            
            temperatures[f"q:_{heat_flux}"][f"alpha_{alpha_this}"] = {}
            
            # iterate over each time step
            for t in t_grid_this:
                
                # create vector b
                b = vector_b(bc_surface, bc_back, upsilon_this, space_divisions, dx, k, T, T_initial, q_array)
                
                # calculate value of future temperature
                Tn = np.linalg.solve(A,b)
            
                # update present temperature
                T = Tn.copy()

                # store temperature profile at this time in the overall dictionary
                temperatures[f"q:_{heat_flux}"][f"alpha_{alpha_this}"][f"t_{t}"] = Tn
            
    return temperatures   
    
    

# function to create tri-diagonal matrix
def tridiag_matrix(bc_surface, bc_back, upsilon, space_divisions, dx, k):
    """
    Creates tridiagonal matrix A
    Linear system to be solved is Ax = b, and x represents temperature values at time n+1

    Parameters:
    ----------
    bc_surface: boundary condition at the surface
        str
    
    bc_back: boundary condition at the back
        str
    
    upsilon: Fourier number divided by 2. Upsilon = alpha*dt/2*dx2
        float
        
    space_divisions: number of nodes in the spatial domain
        int
        
    dx: size of cell in space domain in m
        float
        
    k: thermal conductivity in W/mK
        float
    
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
    if bc_surface == "Linear":
        
        h_total = 45
        pass
    
    elif bc_surface == "Non-linear":
        pass
    

#    elif bc == "Robin":
#        
#        h, _ = bc_data
#        A[0,0] = 1 + 2*upsilon + 2*dx*h/k
#        A[0,1] = - 2 * upsilon
    
    # adjust matrix for the back boundary conditions
    if bc_back == "Insulation":
        pass
    
    elif bc_back == "Aluminium_block":
        pass
#    A[-1, -2] = - 2 * upsilon
#    A[-1, -1] = 1 + 2 * upsilon
#
    return A




def vector_b(bc_surface, bc_back, upsilon, space_divisions, dx, k, T, T_initial, q_array):
    """
    Calculates vector b. Right hand side of linear system of equations

    Parameters:
    ----------
    bc_surface: boundary condition at the surface
        str
    
    bc_back: boundary condition at the back
        str
    
    upsilon: Fourier number divided by 2. Upsilon = alpha*dt/2*dx2
        float
        
    space_divisions: number of nodes in the spatial domain
        int
        
    dx: size of cell in space domain in m
        float
        
    k: thermal conductivity in W/mK
        float
        
    T: array of present temperatures
        np.array

    T_initial: initial temperature in K
        int
        
    q_array: array of size t_grid that contains the incident heat flux at each time step
        np.array
    
    Returns:
    -------
    b: vector to solve linear system of equations
        np.array
    """
    
    # matrix B, similar to matrix A but multiplies T at present
#    B = np.diagflat([upsilon for i in range(space_divisions - 1)], -1) +\
#        np.diagflat([1 - 2 * upsilon for i in range(space_divisions)]) +\
#        np.diagflat([upsilon for i in range(space_divisions - 1)], 1)

    # Calculate vector b
    b = np.zeros(space_divisions)
#    b[1:-1] = B[1:-1, :].dot(T)
#
#    # adjust b vector depending on the surface
#    if bc == "Dirichlet":
#        
#        T_surface = bc_data
#        b[0] = T_surface
#        
#    elif bc == "Neunman":
#        
#        q = bc_data
#        b[0] = 2*upsilon*T[1] + (1 - 2*upsilon)*T[0] + (4*dx*q*upsilon)/(k)
#        
#    elif bc == "Robin":
#        
#        h, T_gas = bc_data
#        b[0] =  2*upsilon*T[1] + (1 - 2*upsilon - 2*dx*h/k)*T[0] + 4*dx*h*T_gas/k
#    
#    # adjust b vector to the back boundary conditions (semi-infinite)
#    b[-2] = T_initial
#    b[-1] = T_initial
    
    return b    
    



