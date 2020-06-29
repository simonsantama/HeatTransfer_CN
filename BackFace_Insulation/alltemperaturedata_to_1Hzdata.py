"""
This script takes the temperature data calculated on different time grids (as delta t depends on the thermal diffusivity)
and it outputs a similar dictionary but with temperature data at 1 Hz for all conditions.
"""

# import libraries
import time
import pickle
import numpy as np

# import data created in main_validation.py
with open('total_data_general_backinsulated.pickle', 'rb') as handle:
    total_data_general = pickle.load(handle)

# extract all the data from the pickle file
Temperatures = total_data_general["Temperatures"]


# for plots and animations, all data must be logged at the same frequency. This loop creates a similar encompassing dictionary with data at 1 Hz
start= time.time()
Temperatures_1Hz = {}
for level1_hftype in Temperatures:
    
    Temperatures_1Hz[level1_hftype] = {}
    for level2_bcsurface in Temperatures[level1_hftype]:
        
        Temperatures_1Hz[level1_hftype][level2_bcsurface] = {}
        for level3_hf in Temperatures[level1_hftype][level2_bcsurface]:
            
            Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf] = {}
            for level4_alpha in Temperatures[level1_hftype][level2_bcsurface][level3_hf]:
                
                Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][level4_alpha] = {}
                # order time stamps (since dictionary keys are not ordered)
                ordered_timestamps = [float(x.split("_")[1]) for x in 
                                      Temperatures[level1_hftype][level2_bcsurface][level3_hf][level4_alpha].keys()]
                ordered_timestamps.sort()

                # extract the index that is closes to the integer time
                time_int = 0
                # takes value of time stamp that is closes (but larger) that the integer time                    
                for t_number, time_stamp in enumerate(ordered_timestamps):
                    
                    # time 0
                    if t_number == 0:
                        Temperatures_1Hz[level1_hftype][level2_bcsurface][
                                level3_hf][level4_alpha][time_int] = Temperatures[level1_hftype][level2_bcsurface][
                                        level3_hf][level4_alpha][f"t_{time_stamp}"]
                    # other times
                    else:    
                        time_int_new = int(time_stamp)
                        if time_int_new > time_int:
                            Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][
                                    level4_alpha][time_int_new] = Temperatures[level1_hftype][level2_bcsurface][
                                        level3_hf][level4_alpha][f"t_{time_stamp}"]
                            time_int = time_int_new
                        else:
                            pass
                        
print(f"Time taken for reducing the temperature data to 1 Hz: {np.round(time.time() - start,2)} seconds")

 
# extract all the data from the pickle file
t_grid = total_data_general["extra_data"]["t_grid"]
x_grid = total_data_general["extra_data"]["x_grid"]
time_total = total_data_general["extra_data"]["time_total"]
alpha = total_data_general["extra_data"]["alpha"]
k = total_data_general["extra_data"]["k"]
q_all = total_data_general["extra_data"]["q_all"]    
    
    
    
# condense all data to be saved including important plotting parameters
total_data_general_1Hz = {"Temperatures": Temperatures_1Hz, "extra_data": {
        "t_grid":t_grid, "x_grid":x_grid, "time_total": time_total, "alpha": alpha, "k": k, "q_all": q_all},
    }
    
with open('temperatures_backinsulated_1Hz.pickle', 'wb') as handle:
    pickle.dump(total_data_general_1Hz, handle)
    print("Plotting data at 1 Hz saved into pickle")
    