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
        surface_temperature, initial_temperature = BC_values
        analytical_sol = surface_temperature + (initial_temperature - surface_temperature) * special.erf(x_grid / (2 * np.sqrt(material["alpha"] * time_duration)))
        title = "ConstantTemp_SemiInf_"

    elif BC_tuple[0] == "const_nhf" and BC_tuple[1] == "semi_inf":
        surface_nhf, initial_temperature = BC_values
        analytical_sol = initial_temperature +\
            (2 * surface_nhf * np.sqrt(material["alpha"] * time_duration / np.pi) / material["k"]) *\
            np.exp((-x_grid ** 2) / (4 * material["alpha"] * time_duration)) -\
            (surface_nhf * x_grid / material["k"]) * special.erfc(x_grid / (2 * np.sqrt(material["alpha"] * time_duration)))
        title = "ConstantNHF_SemiInf_"

    elif BC_tuple[0] == "convection" and BC_tuple[1] == "semi_inf":
        h_convective, initial_temperature, air_temperature = BC_values
        analytical_sol = initial_temperature + (air_temperature - initial_temperature) * (special.erfc(x_grid / (2 * np.sqrt(material["alpha"] * time_duration))) - np.exp((h_convective * x_grid / material["k"]) + (h_convective**2 * material["alpha"] * time_duration / material["k"]**2)) * special.erfc((x_grid / (2 * np.sqrt(material["alpha"] * time_duration))) + (h_convective * np.sqrt(material["alpha"] * time_duration) / material["k"])))
        title = "Convection_SemiInf_"
    else:
        return None

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
    ax.plot(x_grid, T - 273, linewidth=2.5)
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
