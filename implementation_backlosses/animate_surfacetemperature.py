"""
This script creates animations to compare the data.
Animates surface temperature evolution for all the conditions studied.

Must run main.py and alltemperaturedata_to_1Hzdata.py before running this script (requires the file temperatures_backinsulated_1Hz.pickle)

"""

# import libraries
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
import pickle
import numpy as np
import time

# import 1Hz data
with open('temperatures_backinsulated_1Hz.pickle', 'rb') as handle:
    total_data_general_1Hz = pickle.load(handle)

# extract all the data from the pickle file
Temperatures_1Hz = total_data_general_1Hz["Temperatures"]
t_grid = total_data_general_1Hz["extra_data"]["t_grid"]
x_grid = total_data_general_1Hz["extra_data"]["x_grid"]
time_total = total_data_general_1Hz["extra_data"]["time_total"]
alpha = total_data_general_1Hz["extra_data"]["alpha"]
k = total_data_general_1Hz["extra_data"]["k"]
q_all = total_data_general_1Hz["extra_data"]["q_all"]

# plotting parameters
figure_size = (12,6)
figure_title_fs = 15
subplot_title_fs = 13
axis_labels_fs = 12
text_fs = 11
inset_axis_labels_fs = 9
inset_axis_ticks_fs = 8
y_limits = [-80, 800]
y_limits_secondary = [-6,60]
x_limits = [-30,300]
line_styles = ["-", "--", ":", "-."]
line_color = ["firebrick", "royalblue", "seagreen", "black"] # "blueviolet"]
line_width = 1.75
line_width_inset = 1.25
line_color_q = "dimgrey"

Surface_temperatures_1Hz = {}
# create a dictionary with the surface temperature data
for level1_hftype in Temperatures_1Hz:
    Surface_temperatures_1Hz[level1_hftype] = {}
    
    for level2_bcsurface in Temperatures_1Hz[level1_hftype]:
        Surface_temperatures_1Hz[level1_hftype][level2_bcsurface] = {}
        
        for level3_hf in Temperatures_1Hz[level1_hftype][level2_bcsurface]:
            Surface_temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf] = {}
            
            for level4_alpha in Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf]:
                Surface_temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][level4_alpha] = np.zeros(time_total)
                
                # extract surface temperature at every time step
                for time_stamp in range(time_total):
                    Surface_temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][level4_alpha][time_stamp] = \
                        Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][level4_alpha][time_stamp][0]

# animation of surface temperature evolution
for debug,level1_hftype in enumerate(Temperatures_1Hz):
    
    
    # dummy variable to only animate a given type of heat flux
#    if debug in [0,1,2]:
#        continue
    
    hf_type = level1_hftype.split("_")[1]
    print("-------")
    print(f"Creating animations for {hf_type} heat flux")
    print("-------")
    
    start = time.time()
    
    # create plot and format it
    figure_title = f"{hf_type}" + " $\dot{q}''_{inc}$. Insulated back surface\n" + \
        "Linear BC: $\dot{q}''_{net} = \dot{q}''_{inc} - h_T(T_{surf} - T_{\infty})$\n" + \
        "Non-Linear BC: $\dot{q}''_{net} = \dot{q}''_{inc} - h_c(T_{surf} - T_{\infty}) - \epsilon \sigma T_{surf}^4$\n"

    # create figure and format it
    fig, axis = plt.subplots(2,2, sharey = True, sharex = True, constrained_layout = True, figsize = figure_size)
    fig.suptitle(figure_title, fontsize = figure_title_fs)
    for ax in [axis[0,0], axis[1,0]]:
        ax.set_ylabel("$T_{surf}$ [$^\circ$C]", fontsize = axis_labels_fs)
        ax.set_ylim(y_limits)
        ax.set_yticks(np.linspace(0,y_limits[1],5))
    for ax in [axis[1,0], axis[1,1]]:
        ax.set_xlabel("Time [s]", fontsize = axis_labels_fs)
        ax.set_xlim(x_limits)
        ax.set_xticks(np.linspace(0,x_limits[1], 7))
    secondary_axes = []
    for i, ax in enumerate(axis.flatten()):
        ax.grid(color = "gainsboro", linestyle = "--", linewidth = 0.75)
        
        # secondary axis to show heat flux
        ax1 = ax.twinx()
        ax1.set_ylim(y_limits_secondary)
        if i in [0,2]:
            ax1.axes.yaxis.set_ticklabels([])
        else:
            ax1.set_yticks(np.linspace(0,y_limits_secondary[1],5))
            ax1.set_ylabel("Heat Flux [$kW/m^2$]")
        secondary_axes.append(ax1)
        
        
    level2_keys = list(Temperatures_1Hz[level1_hftype].keys())
    for i,level3_hf in enumerate(Temperatures_1Hz[level1_hftype][level2_keys[0]]):
    
        # title of every subplot is the heat flux                
        if hf_type == "Constant":
            hf = level3_hf.split("_")[1].split(".")[0]
            axis.flatten()[i].set_title(f"q$_{i}$ = {hf} kW/m$^2$")
        elif hf_type == "Linear":
            hf = float(level3_hf.split("_")[1])
            axis.flatten()[i].set_title(f"q$_{i}$ = {np.round(hf,2)}$\cdot$t kW/m$^2s$")
        elif hf_type == "Quadratic":
            hf = float(level3_hf.split("_")[1])
            axis.flatten()[i].set_title(f"q$_{i}$ = {'{:.1e}'.format(hf)}$\cdot$t$^2$ kW/m$^2s^2$")
        elif hf_type == "Sinusoidal":
            hf = float(level3_hf.split("_")[1])
            axis.flatten()[i].set_title(f"q$_{i}$ = Sin({'{:.1e}'.format(hf)}$\cdot$t) kW/m$^2s$")      
    
    # add legend to first plot showing the different values of alpha
    custom_lines = []
    for z, alpha_value in enumerate([alpha[0], alpha[-1]]):
        generic_line = Line2D([0],[0], color = line_color[z], lw = line_width, linestyle = line_styles[0])
        custom_lines.append(generic_line)
    legend_alphas = axis[0,0].legend(custom_lines, [alpha[0], alpha[-1]], fancybox = True, 
                        title = r"$\alpha$ [$m^2/s$]", 
                        loc = "upper left", ncol = 1) 
    axis[0,0].add_artist(legend_alphas)

    # add legend to thrid plot showing the different values of k
    custom_lines = []
    for z, alpha_value in enumerate([k[0], k[-1]]):
        generic_line = Line2D([0],[0], color = line_color[z], lw = line_width, linestyle = line_styles[0])
        custom_lines.append(generic_line)
    legend_alphas = axis[0,1].legend(custom_lines, [k[0], k[-1]], fancybox = True, 
                        title = r"$k$ [$W/mK$]", 
                        loc = "upper left", ncol = 1) 
    axis[0,1].add_artist(legend_alphas)  
    
    # add legend to thrid plot showing the different values of k
    custom_lines = []
    for z in range(2):
        generic_line = Line2D([0],[0], color = "black", lw = line_width, linestyle = line_styles[z])
        custom_lines.append(generic_line)
    legend_alphas = axis[1,0].legend(custom_lines, ["Linear", "Non-linear"], fancybox = True, 
                        title = "Surface BC", 
                        loc = "upper left", ncol = 1) 
    axis[1,0].add_artist(legend_alphas) 
    

    # add legend to fourth subplot to show the value of different heat fluxes in the insert (requires dummy lines)
    generic_line = [Line2D([0],[0], color = line_color_q, lw = line_width, linestyle = ":")]
    axis[1,1].legend(generic_line, ["Heat Flux"], 
        fancybox = True, loc = "upper left", ncol = 1)
        
    # add text with counter to the last subplot
    counter = axis[0,0].text(150,670, "Time: 0 seconds", fontsize = text_fs,
                  bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.5'))

    # to the all lines list append the heat flux line plots that are shown in the inset
    t_grid_q = np.arange(0,time_total)
    q_lines = []
    if hf_type == "Constant":
        q_values = q_all[0]
        q_array_plot = [np.zeros_like(t_grid_q) + y for y in q_values]
    elif hf_type == "Linear":
        q_values = q_all[1]
        q_array_plot = [y*t_grid_q for y in q_values]
    elif hf_type == "Quadratic":
        q_values = q_all[2]
        q_array_plot = [y*t_grid_q*t_grid_q for y in q_values]
    elif hf_type == "Sinusoidal":
        amplitude = q_all[3][0]
        mean = q_all[3][1]
        q_values = q_all[3][2]
        q_array_plot = [mean+amplitude*np.sin(y*t_grid_q) for y in q_values]


    # create list to store every plotting line
    all_lines = []
    all_lines_info = []
    level2_keys = list(Temperatures_1Hz[level1_hftype].keys())
    level3_keys = list(Temperatures_1Hz[level1_hftype][level2_keys[0]].keys())
    level4_keys = list(Temperatures_1Hz[level1_hftype][level2_keys[0]][level3_keys[0]].keys())
    
    for j, level3_hf in enumerate(level3_keys):
        for linearity,level2_bcsurface in enumerate(level2_keys):
            for z, alpha_plot in enumerate([level4_keys[0], level4_keys[-1]]):
                line, = axis.flatten()[j].plot([],[], linewidth = line_width, color = line_color[z],
                                    linestyle = line_styles[linearity])
                all_lines.append(line)
                # all_lines_info used to get an idea of what each line corresponds to
                all_lines_info.append(f"{alpha_plot} {level2_bcsurface} {level3_hf} {line_color[z]}")
    
    
    # add lines for plottign the heat flux
    for q_num,q_val in enumerate(q_values):
            q_line, = secondary_axes[q_num].plot([],[], linewidth = line_width_inset, color = line_color_q, 
                                       linestyle = ":")
            all_lines.append(q_line)
            all_lines_info.append(f"{q_val}")
        
               
    # init function for FuncAnimate
    def init():
        for line in all_lines:
                line.set_data([],[])
        counter.set_text("Time: 0 seconds")
        return all_lines
        
    # animate function for FuncAnimate
    def animate(k):
        # update counter
        counter.set_text(f"Time: {k} seconds")
        
        for l,line in enumerate(all_lines):
            
            # up to line 15 is temperatures profiles
            if l < 16:
                
                # lines 0,4,8,12 plot surface linear, lower alpha
                if l%4 == 0:
                    level2_key = "Surface_Linear"
                    level3_key = f"q:_{q_values[int(l/4)]}"
                    level4_key = level4_keys[0]

                    temperature_data = Surface_temperatures_1Hz[level1_hftype][level2_key][level3_key][level4_key]
                    y = temperature_data[:k] - 288
                    line.set_data(t_grid_q[:k],y)
                    
                # lines 1,5,9,13 plot linear surface, higher alpha
                if l%4 == 1:
                    level2_key = "Surface_Linear"
                    level3_key = f"q:_{q_values[int(np.floor(l/4))]}"
                    level4_key = level4_keys[-1]
                    
                    temperature_data = Surface_temperatures_1Hz[level1_hftype][level2_key][level3_key][level4_key]
                    y = temperature_data[:k] - 288
                    line.set_data(t_grid_q[:k],y)
                    
                # lines 2, 6, 10, 14 plot non linear surface, low alpha
                if l%4 == 2:
                    level2_key = "Surface_Non-linear"
                    level3_key = f"q:_{q_values[int(np.floor(l/4))]}"
                    level4_key = level4_keys[0]
                    
                    temperature_data = Surface_temperatures_1Hz[level1_hftype][level2_key][level3_key][level4_key]
                    y = temperature_data[:k] - 288
                    line.set_data(t_grid_q[:k],y)             
                
                # lines 3, 7, 11, 15 plot non linear surface, high alpha
                if l%4 == 3:
                    level2_key = "Surface_Non-linear"
                    level3_key = f"q:_{q_values[int(np.floor(l/4))]}"
                    level4_key = level4_keys[-1]
                    
                    temperature_data = Surface_temperatures_1Hz[level1_hftype][level2_key][level3_key][level4_key]
                    y = temperature_data[:k] - 288
                    line.set_data(t_grid_q[:k],y)
            
            # last four lines are the heat flux plots in the inset
            elif 16<= l:
                x = t_grid_q[:k]
                y = q_array_plot[l-16][:k]
                line.set_data(x,y)
                
        return all_lines

    # create and save animations
    anim = FuncAnimation(fig, animate, init_func=init,
                                   frames=time_total, 
                                   interval=20, blit=True)
    
       
    # save animation in the corresponding folder
    file_name_animation = f"./animations_sufacetemperature/{hf_type}_heatflux-SurfaceTemperature"
    anim.save(f'{file_name_animation}.mp4', dpi = 300, fps = 30)
    plt.close()
        

    print(f" time taken for {hf_type} heat flux animations: {np.round(time.time() - start,2)} seconds")