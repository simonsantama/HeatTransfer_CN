"""
This script creates videos of the comparison of the evolution of the temperature profiles for
the three boundary conditions as calculated in main_validation.py

"""
# import libraries
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle
import numpy as np

# import data 
with open('validation_temperatures.pickle', 'rb') as handle:
    validation_temperatures = pickle.load(handle)

T_analytical = validation_temperatures["analytical"]
T_numerical = validation_temperatures["numerical"]

with open('grid_time.pickle', 'rb') as handle:
    grid_time = pickle.load(handle)
    
x_grid = grid_time[0]
time_total = grid_time[1] - 1
    


# plot title description
bc_description = ["Constant surface temperature = 500 $^\circ$C",
                  "Constant surface heat flux = 10 kW/m2",
                  "Convection h = 10 W/m2K - T_gas = 500 $^\circ$C"]

for i, bc in enumerate(T_analytical):
    
    # create different animations for each boundary condition
    print(f"Animation {bc} boundary condition")
    
    # create figure and set style
    fig, ax = plt.subplots(1,1)
    ax.set_xlabel("Sample length [mm]")
    ax.set_xlim([-10,110])
    ax.set_xticks([0,20,40,60,80,100])
    ax.set_ylabel("Temperature [$^\circ$C]")
    ax.set_ylim([-50,550])
    ax.set_yticks([0,100,200,300,400, 500])
    ax.grid(color = "gainsboro", linestyle = "--", linewidth = 0.75)
    ax.set_title(f"{bc} boundary condition.\n {bc_description[i]}")
    
    # initialize plotted elements
    text = ax.text(65, 350, "time = 0 seconds", fontsize = 12)
    line_analytic, = ax.plot([],[], linewidth = 2, color = "maroon", label = "Analytical solution")
    line_numeric, = ax.plot([],[], linewidth = 0, color = "royalblue", marker = "o", markersize = 3, 
                            markerfacecolor = "lightskyblue", label = "Numerical solution")
    ax.legend(loc = "upper right")
    
    # init function for FuncAnimation
    def init():
        line_analytic.set_data([],[])
        line_numeric.set_data([],[])
        text.set_text("time = 0 seconds")
        return line_analytic, line_numeric, text
    
    # animate function for Funcanimation
    def animate(i):
        x = x_grid*1000
        y_analytic = T_analytical[bc][i+1] - 300
        y_numeric = T_analytical[bc][i+1] - 300
        
        # only update the time counter once every ten seconds
        text.set_text(f"time = {i} seconds")

        line_analytic.set_data(x, y_analytic)
        line_numeric.set_data(x[::10],y_numeric[::10])
        return line_analytic, line_numeric, text
    
    anim = FuncAnimation(fig, animate, init_func=init,
                                    frames = 100,
#                                   frames=time_total, 
                                   interval=20, blit=True)


    anim.save(f'{bc}.mp4', dpi = 300, fps = 30)
    
