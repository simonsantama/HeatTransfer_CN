"""
Crank-Nicholson implementation of 1D heat transfer with simplified boundary conditions for comparison with
analytical solutions.
"""

# import libraries
import numpy as np


def cn_solver(x_grid, t_grid, upsilon, bc, space_divisions):
    """
    Applies Crank-Nicolson method to numerically solve the heat diffusion over the specified domain
    
    Validation function: Exposed surface uses dirichlet, neunman or robin boundary conditions.
    Semi-infinite solid
    
    Parameters:
    ----------
    x_grid: one dimensional spatial domain
        np.array
        
    t_grid: temporal domain
        np.array
        
    upsilon: Fourier number divided by 2
        float
        
    bc: boundary condition
        str
        
    space_divisions: number of nodes in the spatial domain
        int    
    
    Returns:
    -------
    Temperature: dictionary with array of temperatures for each time
        np.array
    
    """
    
    # create dictionary to store all temperatures at all times
    Temperature = {}
    
    # create matrix A
    A = tridiag_matrix(bc, upsilon, space_divisions)
    print(A)
    
    
    
    Temperature = 100000*x_grid
    
    return Temperature


# function to create tri-diagonal matrix
def tridiag_matrix(bc, upsilon, space_divisions):
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
    if bc == "Dirichlet"
        A[0, 0] = 1
        A[0, 1] = 0
        
    # adjust matrix for the boundary condition at the unexposed surface
    A[-1, -2] = 0
    A[-1, -1] = 1
#
#    if BC_tuple[0] == "const_nhf":
#        A[0, 1] = -2 * upsilon
#        if BC_tuple[1] == "semi_inf":
#            A[-1, -2] = 0
#            A[-1, -1] = 1
#
#    if BC_tuple[0] == "convection":
#        h_convective = BC_values[0]
#        A[0, 0] = (
#            1 + 2 * upsilon +
#            dx * upsilon * h_convective / material["k"] +
#            4 * material["emissivity"] * stefan_boltzman * T[0]**3
#        )
#        A[0, 1] = -2 * upsilon
#        if BC_tuple[1] == "semi_inf":
#            A[-1, -2] = 0
#            A[-1, -1] = 1
#

    return A


