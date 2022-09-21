

# Layered Green & Ampt with redistribution (LGAR) in python

**Description**:  LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi arid climates. LGAR's main selling point is that it closely mimics precipitation partitioning results as simulated by the Richards / Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi arid areas. LGAR as implemented in Python is BMI compatible. LGAR is currently being developed in C as well.

Other things to include:

  - **Status**:  As of this writing (21 September 2022), this model will soon be submitted to WRR.
  - **Links to production or demo instances** Animations of LGAR simulations can be found at: https://www.hydroshare.org/resource/46747d77d0ce4995b1e1d58384e15a09/


## Dependencies

This project requires Python 3.

## Usage

In order to successfully run the LGAR model, the config file must be filled out correctly, and the correct files must be present in the parameter_files, forcing_data_files, and config_files folders. These components are described here.

In the config file for a given LGAR model run, "output_file_name_fluxes" is the file path and name of the output file, which will contain fluxes per time step, cumulative fluxes, forcing data, and mass balance calculations.

"forcing_data_folder" is the path that contains the relevant forcing data files for the model run.

"raw_forcing_data_file_name" is the name of the raw, unformatted forcing data file.

"formatted_forcing_data_name" will be the name of the output, correctly resampled and formatted csv file that LGAR saves to the disk. While this is not necessary for the model to run per se, it is useful to have the reformatted forcing data saved.

"params_file" is the path and name of the parameters file to be used with the model.

"forcing_data_formatter_file" is the path and name of the .py file that converts the raw forcing data file to a format that is usable by LGAR.


The model will run correctly if these variable names are filled out in the config file, and if the files occur in the correct places. LGAR also requires the files LGAR_compute.py, test_env.py, and BMI_base_class.py, to be in the same directory. In order to run LGAR, ensure that the correct config file is indicated in test_env.py, and then navigate in a terminal to the directory containing these files and run "python test_env.py".

The Jupyter notebooks (in vis_files) are useful for visualization of results. HYDRUS_files contains HYDRUS-1D model runs which are set up to simulate the same soil hydraulic conditions and forcing data as various LGAR runs.

## How to test the software

In the outputs folder, there are 6 complete simulation outputs, including 3 simulations of USDA SCAN sites, and 3 simulations with synthetically generated forcing datasets. All of the files necessary to run these simulations are also included. In order to check if LGAR is working properly for you, you can run these simulations and compare your results agains the outputs stored in this repo.

## Known issues



## Getting involved




----




----

## Credits and references

LGAR is
