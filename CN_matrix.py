# ------- CREATES THE MATRICES FOR THE CRANK-NICHOLSON SCHEME AND OUTPUTS THE TEMPERATURE AT EVERY TIME STEP -----

import numpy as np
np.set_printoptions(precision=3)

# --- Define tridiagonal matrix A


def tridiag_matrix(sigma, space_mesh_division, BC_tuple, BC_values, material, dx):
    """
    Creates tridiagonal matrix A

    BC_tuple = (BC_surface, BC_back)

    BC_surface: "const_temp", "const_nhf", "convection", "conv_rad"
    BC_back: "semi_inf" or "al_block"
    """

    A = np.diagflat([-sigma for i in range(space_mesh_division - 1)], -1) +\
        np.diagflat([1 + 2 * sigma for i in range(space_mesh_division)]) +\
        np.diagflat([-sigma for i in range(space_mesh_division - 1)], 1)

    # Boundary conditions (plotting of analytical solutions)
    if BC_tuple[0] == "const_temp":
        A[0, 0] = 1
        A[0, 1] = 0
        if BC_tuple[1] == "semi_inf":
            A[-1, -2] = 0
            A[-1, -1] = 1

    if BC_tuple[0] == "const_nhf":
        A[0, 1] = -2 * sigma
        if BC_tuple[1] == "semi_inf":
            A[-1, -2] = 0
            A[-1, -1] = 1

    h_convective = BC_values[0]
    if BC_tuple[0] == "convection":
        A[0, 0] = 1 + 2 * sigma + (2 * dx * h_convective * sigma) / (material["k"])
        A[0, 1] = -2 * sigma
        if BC_tuple[1] == "semi_inf":
            A[-1, -2] = 0
            A[-1, -1] = 1

    # Boundary condition 4: convection + radiation
    if BC_tuple[0] == "conv_rad":
        pass

    return A

# --- Define b_vector


def vector_b(sigma, space_mesh_division, T, BC_tuple, BC_values, dx, material):
    """
    Calculates vector at the right hand side of the algebraic equation

    BC_surface: "const_temp", "const_nhf", "convection", "conv_rad"
    BC_back: "semi_inf", "insulated" "al_block"
    """
    B = np.diagflat([sigma for i in range(space_mesh_division - 1)], -1) +\
        np.diagflat([1 - 2 * sigma for i in range(space_mesh_division)]) +\
        np.diagflat([sigma for i in range(space_mesh_division - 1)], 1)

    # Calculate vector b
    b = np.zeros(space_mesh_division)
    b[1:-1] = B[1:-1, :].dot(T)

    # Boundary conditions 1: constant temperature
    if BC_tuple[0] == "const_temp" and BC_tuple[1] == "semi_inf":
        surface_temperature, initial_temperature = BC_values
        b[0] = surface_temperature
        b[-1] = initial_temperature

    # Boundary conditions 2: constant temperature
    if BC_tuple[0] == "const_nhf" and BC_tuple[1] == "semi_inf":
        surface_nhf, initial_temperature = BC_values
        b[0] = (4 * sigma * dx * surface_nhf) / material["k"] + (1 - 2 * sigma) * T[0] + 2 * sigma * T[1]
        b[-1] = initial_temperature

    # Boundary conditions 3: convection
    if BC_tuple[0] == "convection" and BC_tuple[1] == "semi_inf":
        h_convective, initial_temperature, air_temperature = BC_values
        b[0] = 2 * sigma * T[1] + (1 - 2 * sigma - (2 * dx * h_convective * sigma / material["k"])) * T[0] + (4 * dx * h_convective * sigma / material["k"]) * air_temperature
        b[-1] = initial_temperature

    # Boundary condition 4: convection + radiation

    return b
