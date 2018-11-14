# -------  VERIFIES MODEL BY PLOTTING ANALYTICAL AND NUMERICAL SOLUTIONS  -------

import numpy as np
import matplotlib.pyplot as plt
from scipy import special


def verify_plot(Tn, x_grid, x_lim_other, y_lim_other, time_duration, BC_tuple, BC_values, material):
    """
    Plots the analytical and numerical solution to verify the model at time = t_end
    """
    if BC_tuple[0] == "const_temp" and BC_tuple[1] == "semi_inf":
        analytical_sol = BC_values[0] + (BC_values[1] - BC_values[0]) * special.erf(x_grid / (2 * np.sqrt(material["alpha"] * time_duration)))
        title = "ConstantTemp_SemiInf_"
    elif BC_tuple[0] == "const_nhf" and BC_tuple[1] == "semi_inf":
        pass
    elif BC_tuple[0] == "const_convection" and BC_tuple[1] == "semi_inf":
        pass

    fig = plt.figure()
    plt.scatter(x_grid, Tn, s=80, c="k", marker="o", alpha=0.5, label="Model")
    plt.plot(x_grid, analytical_sol, label="Analytical solution")
    plt.xlim(x_lim_other)
    plt.ylim(y_lim_other)
    plt.legend()
    plt.savefig("Verification_" + title + "_" + str(time_duration) + "s.pdf")
    plt.title(title + str(time_duration) + "s")

    plt.show()

    return None
