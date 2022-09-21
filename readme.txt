In order to successfully run the LGAR model, the config file must be filled out correctly, and the correct files must be present in the parameter_files, forcing_data_files, and config_files folders. These components are described here.

In the config file for a given LGAR model run, "output_file_name_fluxes" is the file path and name of the output file, which will contain fluxes per time step, cumulative fluxes, forcing data, and mass balance calculations. 

"forcing_data_folder" is the path that contains the relevant forcing data files for the model run. 

"raw_forcing_data_file_name" is the name of the raw, unformatted forcing data file. 

"formatted_forcing_data_name" will be the name of the output, correctly resampled and formatted csv file that LGAR saves to the disk. While this is not necessary for the model to run per se, it is useful to have the reformatted forcing data saved. 

"params_file" is the path and name of the parameters file to be used with the model. 

"forcing_data_formatter_file" is the path and name of the .py file that converts the raw forcing data file to a format that is usable by LGAR. 


The model will run correctly if these variable names are filled out in the config file, and if the files occur in the correct places. LGAR also requires the files LGAR_compute.py, test_env.py, and BMI_base_class.py, to be in the same directory. In order to run LGAR, ensure that the correct config file is indicated in test_env, and then navigate in a terminal to the directory containing these files and run "python test_env.py". 

The Jupyter notebooks (in vis_files) are useful for visualization of results. HYDRUS_files contains HYDRUS-1D model runs which are set up to simulate the same soil hydraulic conditions and forcing data as various LGAR runs. 