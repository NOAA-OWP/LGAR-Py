###this is the config file for the LGAR model



#######
### model run options

time_step = 300/3600 #this is the time step in hours. The default value is 300/3600, or 5 minutes expressed in hours.

initial_psi = 20000 #The initial capillary head (mm) set everywhere in the soil. Note that LGAR deals with absolute values of psi, so a value here of 1000 physically corresponds to a pressure head of -1000 mm everywhere.

verbose = True #prints many model state, fluxes, and other interesting quantities if set to True

length_of_simulation = 100000 #This is the length of the simulation, in time steps.

#by default, the LGAR model does not save soil moisture profiles. This is achieved with the following line.
time_steps_to_record_profile = ()

###in the event that you want to record soil moisture profiles for each time step, please either manually enter the time steps you want recorded, or use the np.array code to record all time steps.
# time_steps_to_record_profile = (100,200,300)
# import numpy as np

#######





#######
###file names and paths
###this will be the path and name of the output file which contains the fluxes for each time step
output_file_name_fluxes = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/outputs/output_Bushland.pkl'

###this is the path of the directory that contains the forcing data
forcing_data_folder = '~/desktop/alt-modular/Modules/GAR/LGAR_BMI/forcing_data_files/Bushland/'

###this is the name of the raw, unformatted forcing data file
raw_forcing_data_file_name = 'forcing_data_Bushland_raw.csv'

###this is the name of the resampled, correctly formatted forcing data file that is output by the forcing data reformatter
formatted_forcing_data_name = 'forcing_data_resampled_Bushland.csv'

###this is the path and name of the parameters file
params_file = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/parameter_files/params_Bushland.py'

###this is the path and name of the file the reformats raw forcing data
forcing_data_formatter_file = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/forcing_data_files/Bushland/forcing_data_formatter_Bushland.py'
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
