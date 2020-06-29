"""
Crank-Nicholson implementation of 1D heat transfer with simplified boundary conditions for comparison with
analytical solutions.
"""

# import libraries
import numpy as np


def cn_solver(x_grid, t_grid, upsilon, bc, bc_data, space_divisions, T_initial, dx, k):
    """
    Applies Crank-Nicolson method to numerically solve heat diffusion over the specified domain
    
    Validation function: Exposed surface uses dirichlet, neunman or robin boundary conditions.
    Semi-infinite solid
    
    Parameters:
    ----------
    x_grid: one dimensional spatial domain
        np.array
        
    t_grid: temporal domain
        np.array
        
    upsilon: Fourier number divided by 2 = alpha*dt/2*dx2
        float
        
    bc: boundary condition
        str
        
    bc_data: constant surface temperature (Dirichlet), surface heat flux (Neunman), h and T_gas (Robin)
        int or tuple    
        
    space_divisions: number of nodes in the spatial domain
        int
        
    T_initial: initial temperature in K
        int

    dx = size of cell in space domain in m
        float
        
    k = thermal conductivity in W/mK
        float

    Returns:
    -------
    Temperature: dictionary with array of temperatures for each time
        np.array
    
    """
    
    # create dictionary to store all temperatures at all times
    Temperature = {}
    
    # create matrix A
    A = tridiag_matrix(bc, upsilon, space_divisions, bc_data, dx, k)
    
    # initialise temperature arrays for present and future temperatures
    T = np.zeros_like(x_grid) + T_initial;
    Tn = np.zeros_like(x_grid)
    
    for dt in t_grid:
        
        # create vector b
        b = vector_b(bc, space_divisions, upsilon, T, bc_data, T_initial, dx, k)
        
        # calculate value of future temperature
        Tn = np.linalg.solve(A,b)
        
        # append temperature values to dictionary
        Temperature[dt] = Tn
        
        # update present temperature
        T = Tn.copy()
    
    return Temperature


# function to create tri-diagonal matrix
def tridiag_matrix(bc, upsilon, space_divisions, bc_data, dx, k):
    """
    Creates tridiagonal matrix A
    Linear system to be solved is Ax = b, and x represents temperature values at time n+1

    Parameters:
    ----------
    bc: boundary condition
        str
    
    upsilon: Fourier number divided by 2 = alpha*dt/2*dx2
        float
        
    space_divisions: number of nodes in the spatial domain
        int

    bc_data: constant surface temperature (Dirichlet), surface heat flux (Neunman), h and T_gas (Robin)
        int or tuple 
        
    dx = size of cell in space domain in m
        float
        
    k = thermal conductivity in W/mK
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
    if bc == "Dirichlet":
        A[0, 0] = 1
        A[0, 1] = 0
    elif bc == "Neunman":
        A[0,0] = 1 + 2*upsilon
        A[0,1] = - 2 * upsilon
    elif bc == "Robin":
        h, _ = bc_data
        A[0,0] = 1 + 2*upsilon + 2*dx*h/k
        A[0,1] = - 2 * upsilon
        
    # adjust matrix for the boundary condition at the unexposed surface (insulated)
    A[-1, -2] = - 2 * upsilon
    A[-1, -1] = 1 + 2 * upsilon

    return A


# vector b
def vector_b(bc, space_divisions, upsilon, T, bc_data, T_initial, dx, k):
    """
    Calculates vector b. Right hand side of linear system of equations

    Parameters:
    ----------
    bc: boundary condition
        str

    space_divisions: number of nodes in the spatial domain
        int
        
    upsilon: Fourier number divided by 2 = alpha*dt/2*dx2
        float
        
    T: array of present temperatures
        np.array
        
    bc_data: constant surface temperature (Dirichlet), surface heat flux (Neunman), h and T_gas (Robin)
        int or tuple
        
    T_initial: initial temperature in K
        int
        
    dx = size of cell in space domain in m
        float
        
    k = thermal conductivity in W/mK
        float
    
    Returns:
    -------
    b: vector to solve linear system of equations
        np.array
    """
    
    # matrix B, similar to matrix A but multiplies T at present
    B = np.diagflat([upsilon for i in range(space_divisions - 1)], -1) +\
        np.diagflat([1 - 2 * upsilon for i in range(space_divisions)]) +\
        np.diagflat([upsilon for i in range(space_divisions - 1)], 1)

    # Calculate vector b
    b = np.zeros(space_divisions)
    b[1:-1] = B[1:-1, :].dot(T)

    # adjust b vector depending on the surface
    if bc == "Dirichlet":
        T_surface = bc_data
        b[0] = T_surface
        
    elif bc == "Neunman":
        q = bc_data
        b[0] = 2*upsilon*T[1] + (1 - 2*upsilon)*T[0] + (4*dx*q*upsilon)/(k)
        
    elif bc == "Robin":
        h, T_gas = bc_data
        b[0] =  2*upsilon*T[1] + (1 - 2*upsilon - 2*dx*h/k)*T[0] + 4*dx*h*T_gas/k
    
    # adjust b vector to the back boundary conditions (semi-infinite)
    b[-2] = T_initial
    b[-1] = T_initial
    
    return b    
    



