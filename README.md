

# Layered Green & Ampt with redistribution (LGAR) in python

contact: peter.lafollette@noaa.gov, plafollette@lynker.com

**Description**:  LGAR is a model which partitions precipitation into infiltration and runoff, and is designed for use in arid or semi arid climates. LGAR's main selling point is that it closely mimics precipitation partitioning results as simulated by the Richards / Richardson equation (RRE), without the inherent reliability and stability challenges the RRE poses. Therefore, this model is useful when accurate, stable precipitation partitioning simulations are desired in arid or semi arid areas. LGAR as implemented in Python is BMI compatible. LGAR is currently being developed in C as well.

  - **Status**:  A manuscript describing LGAR was submitted to WRR in October 2022.
  - **Links to production or demo instances** Animations of LGAR simulations can be found at: https://www.hydroshare.org/resource/46747d77d0ce4995b1e1d58384e15a09/


## Dependencies

This project requires Python 3.

## Usage

In order to successfully run the LGAR model, the config file must be filled out correctly, and the correct files must be present in the parameter_files, forcing_data_files, and config_files folders.



## Config file contents

The contents of the first two blocks of the config file for a model run are described here.

First block:

time_step: this is the model time step, expressed in hours. It defaults to a value of 300/3600, or 5 minutes expressed in hours.

initial_psi: this is the uniform capillary head throughout the model domain expressed in mm. Note that LGAR uses absolute values for capillary head, such that a value of 20000 mm for initial_psi physically represents soil with a capillary head of -20000 mm.

verbose: this can be True or False, where no output is printed to the screen during a model run if it is False.

length_of_simulation: This is the length of the simulation, in time steps.

closed_form_capillary_drive_term: set to True or False. This determines if the capillary drive term G is calculated with the numeric integral of hydraulic conductivity with respect to capillary head, or if the equivalence relations between van Genuchten and Brooks - Corey hydraulic models is used to achieve an approximate, closed form equation for G. Setting this value to True generally significantly increases the speed while insignificantly altering model results.  

Second block:

output_file_name_fluxes: this will be the path and name of the output file which contains the fluxes for each time step of the simulation

params_file: this is the path and name of the parameters file

forcing_data_file: this is the forcing data file that is in the correct format for LGAR-Py


## parameter_files

Each parameter file, in the parameter_files folder, only has to be edited in the first block, which contains options related to soil hydraulic parameters, number of layers, maximum ponded head, and options related to the PET -> AET correction.


## forcing_data_files

This folder contains sub folders for each model run. All that is necessary to run LGAR is correctly formatted forcing data as a .csv. Raw USDA SCAN data and notebooks that convert these raw data to the format usable by LGAR are also provided. Currently, it is necessary that the forcing data resolution is the same as time_step as specified in the config file.

## other useful files

LGAR also requires the files LGAR_compute.py, test_env.py, and BMI_base_class.py, to be in the same directory.

The Jupyter notebooks (in vis_files) are useful for visualization of results. HYDRUS_files contains HYDRUS-1D model runs which are set up to simulate the same soil hydraulic conditions and forcing data as various LGAR runs.


## Build instructions
git clone https://github.com/NOAA-OWP/LGAR-Py -b LGAR-Py_public


## Running the model

In order to run LGAR, ensure that the correct config file is indicated in test_env.py, and then navigate in a terminal to the directory containing test_env.py and enter "python test_env.py".



## How to test the software

In the outputs folder, there are 6 complete simulation outputs, including 3 simulations of USDA SCAN sites, and 3 simulations with synthetically generated forcing datasets. All of the files necessary to run these simulations are also included. In order to check if LGAR is working properly for you, you can run these simulations and compare your results agains the outputs stored in this repo.

## Limitations

LGAR is designed for arid and semi arid areas only. Specifically, LGAR is only suitable for environments in which PET>precipitation, and in which fluxes to groundwater are negligibly small compared to the other external fluxes (AET and precipitation). While there is technically a flux by which water leaves the model domain through the lower boundary, this is only based on a choice to ensure mass balance closure; namely, that the portion of a wetting front that exceeds the lower boundary is considered to be lost through the lower boundary. Typically, these fluxes tend to be minuscule, especially considering that wetting fronts that reach the lower boundary of the deepest layer are typically very slow. Essentially, LGAR should not be used for groundwater simulations, and should only be used in semi arid or arid areas. 

## Known issues


## Getting involved




----




----

## Credits and references
