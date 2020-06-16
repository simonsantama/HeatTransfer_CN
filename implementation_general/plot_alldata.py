"""
This script creates animations to compare the data

"""

# import libraries
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pickle
import numpy as np
import time
import pandas as pd

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

# plotting parameters
figure_size = (12,6)
figure_title_fs = 15
subplot_title_fs = 13
axis_labels_fs = 12
text_fs = 11
y_limits = [-50, 500]
x_limits = [- x_grid[-1]*1000/10, x_grid[-1]*1000]
line_styles = ["-", "--", ":", "-."]
line_color = ["firebrick", "royalblue", "seagreen", "black"] # "blueviolet"]
line_width = 1.75

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
                    for time_int in range(time_total):
                        for time_stamp in Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf][level5_alpha]:
                            if float(time_stamp.split("_")[1]) > time_int:
                                Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][level4_hf][level5_alpha][time_int] = Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf][level5_alpha][time_stamp]

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
                      
            # create list to store every plotting line
            all_lines = []
            for i,level4_hf in enumerate(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback]):
                
                # title of every subplot is the heat flux                
                if hf_type == "Constant":
                    hf = level4_hf.split("_")[1].split(".")[0]
                    axis.flatten()[i].set_title(f"q = {hf} kW/m$^2$")
                elif hf_type == "Linear":
                    hf = int(level4_hf.split("_")[1].split(".")[0])/100
                    axis.flatten()[i].set_title(f"q = {hf}*t kW/m$^2s$")

                # initialize plotted elements
                for j, level5_alpha in enumerate(Temperatures[level1_hftype][level2_bcsurface][level3_bcback][level4_hf]):
            
                        line, = axis.flatten()[i].plot([],[], linewidth = line_width, color = line_color[j], 
                                               linestyle = line_styles[j], label = "{:.1e}".format(alpha[j]))
                        all_lines.append(line)
                        
                
                # collect all the data into the temperature_data_animation numpy array
                
            # add text with counter to first subplot and legend to second    
            axis[0,1].legend(fancybox = True, title = r"$\alpha$", loc = "upper right")
            counter = axis[0,0].text(15,400, "time: 0 seconds", fontsize = text_fs)
            
            # dictionary keys used to access the temperature data with the animate function
            level3_keys = list(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback].keys())
            
            
            # init function for FuncAnimate
            def init():
                for line in all_lines:
                        line.set_data([],[])
                counter.set_text("time: 0 seconds")
                return all_lines
            
            # animate function for FuncAnimate
            def animate(k):
                # update counter
                counter.set_text(f"time: {k} seconds")
                # plot in mm
                x = x_grid*1000
                
                level4_keys = list(Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][level3_keys[0]].keys())
                for l,line in enumerate(all_lines):
                    
                    # first four lines are the first subplot
                    if l < 4:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[0]][level4_keys[l]]
                        
                        # plot in C. data is in K.
                        y = temperature_data[k]
                        line.set_data(x,y)
                    
                    # second subplot
                    elif 4 <= l < 8:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[0]][level4_keys[l-4]]
                        y = temperature_data[k]
                        line.set_data(x,y)
                        
                    # third subplot
                    elif 8 <= l < 12:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[0]][level4_keys[l-8]]
                        y = temperature_data[k]
                        line.set_data(x,y)
                    
                    # fourth subplot
                    elif 12 <= l:
                        temperature_data = Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_bcback][
                                level3_keys[0]][level4_keys[l-12]]
                        y = temperature_data[k]
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
    

