###this is the config file for the LGAR model

###In the config file, the user has to make the following two blocks contain the correct information. The first block, model run options, contains options related to the simulation. the second block, file names and paths, has to contain the desired input and output file names and paths. The third and final block does not have to be changed.

#######
### model run options

time_step = 300/3600 #this is the time step in hours. The default value is 300/3600, or 5 minutes expressed in hours.

initial_psi = 1000 #The initial capillary head (mm) set everywhere in the soil. Note that LGAR deals with absolute values of psi, so a value here of 1000 physically corresponds to a pressure head of -1000 mm everywhere.

verbose = True #prints many model state, fluxes, and other interesting quantities if set to True

length_of_simulation = 144 #This is the length of the simulation, in time steps. Default for Panama is 68910.

#by default, the LGAR model does not save soil moisture profiles. This is achieved with the following line.
time_steps_to_record_profile = ()

###in the event that you want to record soil moisture profiles for each time step, please either manually enter the time steps you want recorded, or use the np.array code to record all time steps.
# time_steps_to_record_profile = (100,200,300)
# import numpy as np
# time_steps_to_record_profile = np.array(range(0,length_of_simulation+1))



#######
###file names and paths
###this will be the path and name of the output file which contains the fluxes for each time step
output_file_name_fluxes = '/Users/peterlafollette/Desktop/LGAR-Py/outputs/output_synth_3.pkl'

###this is the path of the directory that contains the forcing data
forcing_data_folder = '~/desktop/LGAR-Py/forcing_data_files/synth_3/'

###this is the name of the raw, unformatted forcing data file
raw_forcing_data_file_name = 'forcing_data_synth_3_raw.csv'

###this is the name of the resampled, correctly formatted forcing data file that is output by the forcing data reformatter
formatted_forcing_data_name = 'forcing_data_resampled_synth_3.csv'

###this is the path and name of the parameters file
params_file = '/Users/peterlafollette/Desktop/LGAR-Py/parameter_files/params_synth_3.py'

###this is the path and name of the file the reformats raw forcing data
forcing_data_formatter_file = '/Users/peterlafollette/Desktop/LGAR-Py/forcing_data_files/synth_3/forcing_data_formatter_synth_3.py'
#######




#######
###code that loads forcing data and params
import importlib.util
forcing_data_formatter = importlib.util.spec_from_file_location("forcing_data_formatter", forcing_data_formatter_file)
foo = importlib.util.module_from_spec(forcing_data_formatter)
forcing_data_formatter.loader.exec_module(foo)
forcing_data_formatter = foo
forcing_data_file = forcing_data_formatter
forcing_data = forcing_data_formatter.forcing_data_formatter_fxn(path_string=forcing_data_folder, raw_forcing_data_file_name=raw_forcing_data_file_name, formatted_forcing_data_name=formatted_forcing_data_name,freq=time_step)

params = importlib.util.spec_from_file_location("params", params_file)
foo = importlib.util.module_from_spec(params)
params.loader.exec_module(foo)
params = foo
#######
