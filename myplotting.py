# -------  PLOTTING FUNCTIONS FOR THE HEAT TRANSFER MODEL  -------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import special

mpl.rc('font', size=18)
mpl.rc('font', family='Arial')
figure_size = (16, 12)

# Verify the model by plotting against analytical solution (if it exists)


def verify_plot(Tn, x_grid, x_lim_other, y_lim_other, time_duration, BC_tuple, BC_values, material):
    """
    Plots the analytical and numerical solution to verify the model at time = t_end
    """
    if BC_tuple[0] == "const_temp" and BC_tuple[1] == "semi_inf":
        # BC_values[0] = surface temperature and BC_values[1] = intial temperature
        analytical_sol = BC_values[0] + (BC_values[1] - BC_values[0]) * special.erf(x_grid / (2 * np.sqrt(material["alpha"] * time_duration)))
        title = "ConstantTemp_SemiInf_"

    elif BC_tuple[0] == "const_nhf" and BC_tuple[1] == "semi_inf":
        # BC_values[0] = surface nhf and BC_values[1] = initial temperature
        analytical_sol = BC_values[1] +\
            (2 * BC_values[0] * np.sqrt(material["alpha"] * time_duration / np.pi) / material["k"]) *\
            np.exp((-x_grid ** 2) / (4 * material["alpha"] * time_duration)) -\
            (BC_values[0] * x_grid / material["k"]) * special.erf(x_grid / (2 * np.sqrt(material["alpha"] * time_duration)))
        title = "ConstantNHF_SemiInf_"

    elif BC_tuple[0] == "convection" and BC_tuple[1] == "semi_inf":
        pass

    fig = plt.figure()
    plt.scatter(x_grid, Tn, s=80, c="k", marker="o", alpha=0.5, label="C-N Model")
    plt.plot(x_grid, analytical_sol, label="Analytical solution")
    plt.xlim(x_lim_other)
    plt.ylim(y_lim_other)
    plt.legend(loc="upper left")
    plt.savefig("Verification_" + title + "_" + str(time_duration) + "s.pdf")
    plt.title(title + str(time_duration) + "s")
    plt.show()

    return None


def plot_tempgrad(T, x_grid, figure_size, x_lim, y_lim, title, save_format, show):
    """
    Plots temperature gradient at a given time

    Save format = None for no saving
    """
    fig, ax = plt.subplots(figsize=figure_size)
    ax.plot(x_grid, T, linewidth=2.5)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_xlabel("Depth [m]")
    ax.set_ylabel("Temperature [C]")
    ax.set_title(title)

    # Decide whether or not to save the graph
    if save_format == None:
        pass
    else:
        title = title.replace(" ", "")
        fig.savefig(title + "." + save_format, dpi=300)

    # Decide whether or not to show the graph
    if show == "show":
        plt.show()
    elif show == "not-show":
        pass

    return None
