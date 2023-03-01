###this is the config file for the LGAR model

###In the config file, the user has to make the following two blocks contain the correct information. The first block, model run options, contains options related to the simulation. the second block, file names and paths, has to contain the desired input and output file names and paths. The third and final block does not have to be changed.

#######
### model run options

time_step = 300/3600 #this is the time step in hours. The default value is 300/3600, or 5 minutes expressed in hours.

initial_psi = 1000 #100 #The initial capillary head (mm) set everywhere in the soil. Note that LGAR deals with absolute values of psi, so a value here of 1000 physically corresponds to a pressure head of -1000 mm everywhere.

verbose = True #prints many model state, fluxes, and other interesting quantities if set to True

length_of_simulation = 144#12#67#144 #This is the length of the simulation, in time steps.

closed_form_capillary_drive_term = False

#by default, the LGAR model does not save soil moisture profiles. This is achieved with the following line.
time_steps_to_record_profile = ()

###in the event that you want to record soil moisture profiles for each time step, please either manually enter the time steps you want recorded, or use the np.array code to record all time steps.
# time_steps_to_record_profile = (100,200,300)
# import numpy as np
# time_steps_to_record_profile = np.array(range(0,length_of_simulation+1))
#######





#######
###file names and paths
###all the user has to do is specify the relative paths and names of the desired output file that will be created, and the parameters and forcing data files that should exist at the start of the model run
import os, sys
output_path_and_name = 'outputs/output_synth_1.pkl'
params_file_path_and_name = 'parameter_files/params_synth_1.py'
forcing_data_file_path_and_name = 'forcing_data_files/synth_1/forcing_data_resampled_synth_1.csv'

###below this line the user doesn't have to edit anything
output_file_name_fluxes = os.path.join(os.getcwd(),output_path_and_name)
params_file = os.path.join(os.getcwd(),params_file_path_and_name)
forcing_data_file = os.path.join(os.getcwd(),forcing_data_file_path_and_name)
#######







#######
###code that loads forcing data and params
import pandas as pd
import importlib.util
forcing_data = pd.read_csv(forcing_data_file,index_col=0,parse_dates=True)

params = importlib.util.spec_from_file_location("params", params_file)
foo = importlib.util.module_from_spec(params)
params.loader.exec_module(foo)
params = foo
#######
