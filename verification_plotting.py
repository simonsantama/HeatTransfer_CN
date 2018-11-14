import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Plots the analytical and numerical solutions to verify the behaviour of the code (Boundary conditions 1, 2 and 3 as defined in Incropera)

x = sp.symbols("x")  # depth
u = sp.symbols("u")  # temperature


def verify_plot(Tn, x_lim_other, y_lim_other, time_duration, BC_tuple):
    """
    Plots the analytical and numerical solution to verify the model
    """
    if BC_tuple[0] == "const_temp" and BC_tuple[1] == "semi_inf":
        print("so far so good")
    return None
