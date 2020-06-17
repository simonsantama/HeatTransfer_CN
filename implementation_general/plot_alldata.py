"""
This script creates animations to compare the data

"""

# import libraries
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pickle
import numpy as np
import time

# import data created in main_validation.py
with open('total_data_general.pickle', 'rb') as handle:
    total_data_general = pickle.load(handle)

# extract all the data from the pickle file
Temperatures = total_data_general["Temperatures"]
t_grid = total_data_general["extra_data"]["t_grid"]
x_grid = total_data_general["extra_data"]["x_grid"]
time_total = total_data_general["extra_data"]["time_total"]
alpha = total_data_general["extra_data"]["alpha"]
k = total_data_general["extra_data"]["k"]
q_all = total_data_general["extra_data"]["q_all"]

# plotting parameters
figure_size = (12,6)
figure_title_fs = 15
subplot_title_fs = 13
axis_labels_fs = 12
text_fs = 11
inset_axis_labels_fs = 9
inset_axis_ticks_fs = 8
y_limits = [-50, 500]
x_limits = [- x_grid[-1]*1000/10, x_grid[-1]*1000]
line_styles = ["-", "--", ":", "-."]
line_color = ["firebrick", "royalblue", "seagreen", "black"] # "blueviolet"]
line_width = 1.75
line_width_inset = 1.25
line_color_q = "dimgrey"

# for plots and animations, all data must be logged at the same frequency. This loop creates a similar encompassing dictionary with data at 1 Hz
start= time.time()
Temperatures_1Hz = {}
for level1_hftype in Temperatures:
    
    Temperatures_1Hz[level1_hftype] = {}
    for level2_bcsurface in Temperatures[level1_hftype]:
        
        Temperatures_1Hz[level1_hftype][level2_bcsurface] = {}
        for level3_bcback in Temperatures[level1_hftype][level2_bcsurface]:
            
            Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback] = {}
            for level4_hf in Temperatures[level1_hftype][level2_bcsurface][level3_bcback]:
                
                Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][level4_hf] = {}
                for level5_alpha in Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf]:
                    
                    Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][level4_hf][level5_alpha] = {}
                    # order time stamps (since dictionary keys are not ordered)
                    ordered_timestamps = [float(x.split("_")[1]) for x in 
                                          Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf][level5_alpha].keys()]
                    ordered_timestamps.sort()

                    # extract the index that is closes to the integer time
                    time_int = 0
                    # takes value of time stamp that is closes (but larger) that the integer time                    
                    for t_number, time_stamp in enumerate(ordered_timestamps):
                        
                        # time 0
                        if t_number == 0:
                            Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                    level4_hf][level5_alpha][time_int] = Temperatures[level1_hftype][level2_bcsurface][
                                            level3_bcback][level4_hf][level5_alpha][f"t_{time_stamp}"]
                        # other times
                        else:    
                            time_int_new = int(time_stamp)
                            if time_int_new > time_int:
                                Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                    level4_hf][level5_alpha][time_int_new] = Temperatures[level1_hftype][level2_bcsurface][
                                            level3_bcback][level4_hf][level5_alpha][f"t_{time_stamp}"]
                                time_int = time_int_new
                            else:
                                pass
                        
print(f"Time taken for reducing the temperature data to 1 Hz: {np.round(time.time() - start,2)}")


# animation structure that follows the structure of the data
for level1_hftype in Temperatures_1Hz:
    
    hf_type = level1_hftype.split("_")[1]
    print("-------")
    print(f"Creating animations for {hf_type} heat flux")
    print("-------")
    
    start = time.time()
    for level2_bcsurface in Temperatures_1Hz[level1_hftype]:
        bc_surface = level2_bcsurface.split("_")[1]
        
        print(f" surface bc: {bc_surface}")
        
        for level3_bcback in Temperatures_1Hz[level1_hftype][level2_bcsurface]:
            bc_back = level3_bcback.split("_")[1]
            
            print(f"  back bc: {bc_back}")

            figure_title = f"{hf_type} Heat Flux\n Surface BC: {bc_surface}. Back BC: {bc_back}"
            
            # create figure and format it
            fig, axis = plt.subplots(2,2, sharey = True, sharex = True, constrained_layout = True, figsize = figure_size)
            fig.suptitle(figure_title, fontsize = figure_title_fs)
            for ax in [axis[0,0], axis[1,0]]:
                ax.set_ylabel("Temperature [$^\circ$C]", fontsize = axis_labels_fs)
                ax.set_ylim(y_limits)
            for ax in [axis[1,0], axis[1,1]]:
                ax.set_xlabel("Sample depth [mm]", fontsize = axis_labels_fs)
                ax.set_xlim(x_limits)
            for ax in axis.flatten():
                ax.grid(color = "gainsboro", linestyle = "--", linewidth = 0.75)
                
            # add inset to first subplot to show evolution of the heat flux
            ax_inset = inset_axes(axis[0,0], width = "35%", height = "45%", borderpad = 1.5)
            ax_inset.set_xlabel("Time [s]", fontsize = inset_axis_labels_fs)
            ax_inset.set_xlim([0,100])
            ax_inset.set_xticks(np.linspace(0,100,5))
            ax_inset.set_ylabel("q [kW/m$^2$]", fontsize = inset_axis_labels_fs)
            ax_inset.set_ylim([0,60])
            ax_inset.set_yticks(np.linspace(0,60,5))
            ax_inset.grid(color = "gainsboro", linestyle = "--", linewidth = 0.5)
            
            # create list to store every plotting line
            all_lines = []
            for i,level4_hf in enumerate(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback]):
                
                # title of every subplot is the heat flux                
                if hf_type == "Constant":
                    hf = level4_hf.split("_")[1].split(".")[0]
                    axis.flatten()[i].set_title(f"q$_{i}$ = {hf} kW/m$^2$")
                elif hf_type == "Linear":
                    hf = float(level4_hf.split("_")[1])
                    axis.flatten()[i].set_title(f"q$_{i}$ = {hf}$\cdot$t kW/m$^2s$")
                elif hf_type == "Quadratic":
                    hf = float(level4_hf.split("_")[1])
                    axis.flatten()[i].set_title(f"q$_{i}$ = {'{:.1e}'.format(hf)}$\cdot$t$^2$ kW/m$^2s^2$")                    

                # initialize plotted elements
                for j, level5_alpha in enumerate(Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf]):
            
                        line, = axis.flatten()[i].plot([],[], linewidth = line_width, color = line_color[j], 
                                               linestyle = line_styles[j], label = "{:.1e}".format(alpha[j]))
                        all_lines.append(line)
                        
                            
            # add legend to second subplot  
            axis[0,1].legend(fancybox = True, title = r"$\alpha$ [$m^2/s$]", loc = "upper right", ncol = 1)
            # add legend to third subplot but showing the values of k rather than alpha (requires dummy lines)
            custom_lines = []
            for z, k_value in enumerate(k):
                generic_line = Line2D([0],[0], color = line_color[z], lw = line_width, linestyle = line_styles[z])
                custom_lines.append(generic_line)
            axis[1,0].legend(custom_lines, k, fancybox = True, title = "k [$W/mK$]", loc = "upper right", ncol = 1)            
            
            # add legend to fourth subplot to show the value of different heat fluxes in the insert (requires dummy lines)
            custom_lines_q = []
            for z, k_value in enumerate(q_all[0]):
                generic_line = Line2D([0],[0], color = line_color_q, lw = line_width, linestyle = line_styles[z])
                custom_lines_q.append(generic_line)
            axis[1,1].legend(custom_lines_q, ["q$_0$", "q$_1$", "q$_2$", "q$_3$",], 
                fancybox = True, title = "Inset: HF", loc = "upper right", ncol = 1)
            
            # add text with counter to the last subplot
            counter = axis[1,0].text(12,420, "Time: 0 seconds", fontsize = text_fs,
                          bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=0.5'))
            
            # to the all lines list append the heat flux line plots that are shown in the inset
            t_grid_q = np.arange(0,time_total)
            q_lines = []
            if hf_type == "Constant":
                q_values = q_all[0]
                q_array_plot = [np.zeros_like(t_grid_q) + y for y in q_values]
            if hf_type == "Linear":
                q_values = q_all[1]
                q_array_plot = [y*t_grid_q for y in q_values]
            if hf_type == "Quadratic":
                q_values = q_all[2]
                q_array_plot = [y*t_grid_q*t_grid_q for y in q_values]
            
            for z,q_val in enumerate(q_values):
                    q_line, = ax_inset.plot([],[], linewidth = line_width_inset, color = line_color_q, 
                                               linestyle = line_styles[z], label = f"q$_{z}$")
                    all_lines.append(q_line)
            
            
            # dictionary keys used to access the temperature data with the animate function
            level3_keys = list(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback].keys())
            
            
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
                # plot in mm
                x = x_grid*1000
                
                level4_keys = list(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][level3_keys[0]].keys())
                for l,line in enumerate(all_lines):
                    
                    # first four lines are the first subplot
                    if l < 4:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[0]][level4_keys[l]]
                        
                        # plot in C. data is in K.
                        y = temperature_data[k] - 288
                        line.set_data(x,y)
                    
                    # second subplot
                    elif 4 <= l < 8:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[1]][level4_keys[l-4]]
                        y = temperature_data[k] - 288
                        line.set_data(x,y)
                        
                    # third subplot
                    elif 8 <= l < 12:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[2]][level4_keys[l-8]]
                        y = temperature_data[k] - 288
                        line.set_data(x,y)
                    
                    # fourth subplot
                    elif 12 <= l < 16:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[3]][level4_keys[l-12]]
                        y = temperature_data[k] - 288
                        line.set_data(x,y)
                        
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
            file_name_animation = f"./animations_alldata/{hf_type}_heatflux/Surface{bc_surface}_Back{bc_back}"
            anim.save(f'{file_name_animation}.mp4', dpi = 300, fps = 30)
            plt.close()

    print(f"Time taken for {hf_type} heat flux animations: {np.round(time.time() - start,2)} seconds")
