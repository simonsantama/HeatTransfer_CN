"""
This script creates videos to compare the evolution of the temperature profiles for
the three boundary conditions as calculated in main_validation.py

"""

# import libraries
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle
import numpy as np

# import data created in main_validation.py
with open('validation_temperatures.pickle', 'rb') as handle:
    validation_temperatures = pickle.load(handle)

T_analytical = validation_temperatures["analytical"]
T_numerical = validation_temperatures["numerical"]

with open('grid_time.pickle', 'rb') as handle:
    grid_time = pickle.load(handle)
    
x_grid = grid_time[0]
time_total = grid_time[1] - 1
    

# strings to use as ax.plot titles
bc_description = ["Constant surface temperature = 500 $^\circ$C",
                  "Constant surface heat flux = 5 kW/m2",
                  "Convection: h = 50 W/m2K / T_gas = 500 $^\circ$C"]

# create an animation per boundary condition
for i, bc in enumerate(T_analytical):
    
#    if bc == "Dirichlet":
#        continue
    
    # print current boundary condition
    print(f"Animation {bc} boundary condition")
    
    # create figure and set style
    fig, ax = plt.subplots(1,1)
    ax.set_xlabel("Sample depth [mm]")
    
    x_limits = [-20,200]
    y_limits = [[-100,1000],[-100,1000],[-100,1000]]
    ax.set_xlim(x_limits)
    ax.set_xticks(np.linspace(0,x_limits[1] ,6))
    ax.set_ylabel("Temperature [$^\circ$C]")
    ax.set_ylim(y_limits[i])
    ax.set_yticks(np.linspace(0,y_limits[i][1],6))
    
    ax.grid(color = "gainsboro", linestyle = "--", linewidth = 0.75)
    ax.set_title(f"{bc} boundary condition.\n {bc_description[i]}")
    
    # initialize plotted elements
    text = ax.text(115, 650, "time = 0 seconds", fontsize = 12)
    line_analytic, = ax.plot([],[], linewidth = 1.75, color = "maroon", label = "Analytical solution")
    line_numeric, = ax.plot([],[], linewidth = 0, color = "royalblue", marker = "o", markersize = 3.5, 
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
        
        # plot in mm (data is in meters)
        x = x_grid*1000
        
        # plot temperature in C (data is in K)
        y_analytic = T_analytical[bc][i+1] - 300
        line_analytic.set_data(x, y_analytic)
        
        # for numerical results, since we are using the same grid, only plot one out of every 10 points
        y_numeric = T_numerical[bc][i+1] - 300
        line_numeric.set_data(x[::10],y_numeric[::10])
        
        text.set_text(f"time = {i} seconds")

        return line_analytic, line_numeric, text
    
    anim = FuncAnimation(fig, animate, init_func=init,
#                                    frames = 100,
                                   frames=time_total, 
                                   interval=20, blit=True)


    anim.save(f'{bc}.mp4', dpi = 300, fps = 30)
    
