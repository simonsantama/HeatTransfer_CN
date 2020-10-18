"""
This function takes the temperature data calculated on different time grids (as delta t depends on the thermal diffusivity)
and it outputs a similar dictionary but with temperature data at 1 Hz for all conditions.
"""

def data_to1Hz(total_data_general):
    
    # extract all the data from the pickle file
    Temperatures = total_data_general["Temperatures"]
    t_grid = total_data_general["extra_data"]["t_grid"]
    x_grid = total_data_general["extra_data"]["x_grid"]
    time_total = total_data_general["extra_data"]["time_total"]
    alpha = total_data_general["extra_data"]["alpha"]
    k = total_data_general["extra_data"]["k"]
    q_all = total_data_general["extra_data"]["q_all"]     
    
    # extract data at an approximately uniform (1 Hz) logging frequency
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

                    # takes value of time stamp that is closest (but larger) that the integer time                    
                    time_counter = 0
                    for time_stamp in ordered_timestamps:
                        
                        if time_stamp > time_counter:
                            Temperatures_1Hz[level1_hftype][level2_bcsurface][level3_hf][
                                    level4_alpha][time_counter] = Temperatures[level1_hftype][level2_bcsurface][
                                            level3_hf][level4_alpha][f"t_{time_stamp}"]
                            time_counter += 1
                        

    # condense all data to be saved including important plotting parameters
    total_data_1Hz = {"Temperatures": Temperatures_1Hz, 
                      "extra_data": {"t_grid":t_grid, 
                                     "x_grid":x_grid, 
                                     "time_total": time_total, 
                                     "alpha": alpha, 
                                     "k": k, 
                                     "q_all": q_all}}
    
    return total_data_1Hz

    