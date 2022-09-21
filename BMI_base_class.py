###the way I'm choosing to implement BMI with LGAR is to write a base class with all of the BMI functions, which the class LGAR will inherit.
###this way, all of the BMI functions are cleanly separated from LGAR, and this BMI base class might be useful to me in the future.

from datetime import datetime
import numpy as np
from collections import deque
import importlib
import sys
import pandas as pd

class base_BMI():

    def __init__(self):

        ###not sure if all of these are going to be used
        self._values = {}
        self._var_loc = "node"
        self._var_grid_id = 0


        #----------------------------------------------
        # Required, static attributes of the model
        #----------------------------------------------
        self._att_map = {
            'model_name':         'Layered Green-Ampt with redistribution (LGAR)',
            'version':            '1.0',
            'author_name':        'Peter Tazwell La Follette',
            'grid_type':          'scalar',
            'time_step_size':      300/3600,
            'time_units':         'hours' }

        #---------------------------------------------
        # Input variable names (CSDMS standard names)
        #---------------------------------------------
        self._input_var_names = [###names in capital letters are "fake" CSDMS names: I couldn't find CSDMS (CSN) names that fit these variables
            'atmosphere_water__precipitation_volume_flux',
            'land_surface_water__potential_evapotranspiration_volume_flux'] ###note that there is, as far as I can tell, no CSDMS name for evapotranspiration, and only 1 name for transpiration that does not fit this model. Further, the names for evaporation apply to land_surface_water, and this model uses PET from soil water.

        #---------------------------------------------
        # Output variable names (CSDMS standard names)
        #---------------------------------------------
        self._output_var_names = ['land_surface_water__runoff_volume_flux',
                                  'soil_water__infiltration_volume_flux',
                                  'soil_water_sat-zone_top__recharge_volume_flux', #note that this is the lower boundary condition, which is usually 0, but sometimes goes above 0 for individual time steps where part of a wetting front goes deeper than the lower boundary of the deepest layer
                                  'land_surface_water__evapotranspiration_volume_flux', ###note that there is, as far as I can tell, no CSDMS name for evapotranspiration, and only 1 name for transpiration that does not fit this model. Further, the names for evaporation apply to land_surface_water, and this model simulates actual ET from soil water.
                                  'land_surface_water__infiltration_ponding_depth',
                                  'soil_water_wetting-front__depth',
                                  'soil_water__volume_fraction']

        self._var_types = {
                                'atmosphere_water__precipitation_volume_flux':'float32',
                                'land_surface_water__potential_evapotranspiration_volume_flux':'float32',
                                'land_surface_water__runoff_volume_flux':'float32',
                                'soil_water__infiltration_volume_flux':'float32',
                                'soil_water_sat-zone_top__recharge_volume_flux':'float32',
                                'land_surface_water__evapotranspiration_volume_flux':'float32',
                                'land_surface_water__infiltration_ponding_depth':'float32',
                                'soil_water_wetting-front__depth':'float32', #this will be a vector containing at least 3 float64 numbers
                                'soil_water__volume_fraction':'float32' #this will be a vector containing at least 3 float64 numbers
                          }


        #------------------------------------------------------
        # Create a Python dictionary that maps CSDMS Standard
        # Names to the model's internal variable names.
        # This is going to get long,
        #     since the input variable names could come from any forcing...
        #------------------------------------------------------
        self._var_name_units_map = {
                                'atmosphere_water__precipitation_volume_flux':['precip_mm_per_h','mm h-1'],
                                'land_surface_water__potential_evapotranspiration_volume_flux':['PET_mm_per_h','mm h-1'],
                                'land_surface_water__runoff_volume_flux':['runoff_mm_per_h','mm h-1'],
                                'soil_water__infiltration_volume_flux':['actual_infiltration_mm_per_h','mm h-1'],
                                'soil_water_sat-zone_top__recharge_volume_flux':['bottom_flux_mm_per_h','mm h-1'],
                                'land_surface_water__evapotranspiration_volume_flux':['actual_ET_mm_per_h','mm h-1'],
                                'land_surface_water__infiltration_ponding_depth':['ponding_depth_mm ','mm'],
                                'soil_water_wetting-front__depth':['wetting_front_depths','mm'],
                                'soil_water__volume_fraction':['wetting_front_moistures','unitless']
                          }



            # self.precip_mm_per_h = precip_data[i]
            # self.PET_mm_per_h = PET_data[i]
            # self.runoff_mm_per_h = runoff/time_step
            # self.actual_infiltration_mm_per_h = precip_mass_to_add/time_step
            # self.bottom_flux_mm_per_h = free_drainage_demand/time_step
            # self.actual_ET_mm_per_h = actual_ET_demand/time_step
            # self.ponding_depth_mm = self.h_0
            # self.wetting_front_depths #already defined
            # self.wetting_front_moistures #already defined
        #------------------------------------------------------------


    def initialize(self,config_string): #config_string

        #this function takes a string config_string, which is the name of the config file, without the .py on the end. It has to be in the same directory as this .py file to work.


        import importlib.util
        config = importlib.util.spec_from_file_location("config", config_string)
        foo = importlib.util.module_from_spec(config)
        config.loader.exec_module(foo)
        config = foo

        self.mass_balance_error_value_cumulative = 0

        self.precip_mm_per_h = 0
        self.PET_mm_per_h = 0
        self.runoff_mm_per_h = 0
        self.actual_infiltration_mm_per_h = 0
        self.bottom_flux_mm_per_h = 0
        self.actual_ET_mm_per_h = 0
        self.ponding_depth_mm = 0
        self.wetting_front_depths = 0
        self.wetting_front_moistures = np.array([0])

        # self._long_names_to_model_variables_map = {###names here are old
        #                         'atmosphere_water__precipitation_mass_flux':self.precip_mm_per_h,
        #                         'SOIL_WATER__POTENTIAL_EVAPOTRANSPIRATION_MASS_FLUX':self.PET_mm_per_h,
        #                         'land_surface_water__runoff_mass_flux':self.runoff_mm_per_h,
        #                         'soil_water__infiltration_mass_flux':self.actual_infiltration_mm_per_h,
        #                         'soil_water_sat-zone_top__recharge_mass_flux':self.bottom_flux_mm_per_h,
        #                         'SOIL_WATER__ACTUAL_EVAPOTRANSPIRATION_MASS_FLUX':self.actual_ET_mm_per_h,
        #                         'land_surface_water__infiltration_ponding_depth':self.ponding_depth_mm,
        #                         'soil_water_wetting-front__depth':self.wetting_front_depths,
        #                         'soil_water__volume_fraction':self.wetting_front_moistures
        #                   }



        self.current_time = 0
        self.current_time_step = 0
        ####config = importlib.import_module(config_string)
        self.config = config
        #self.folder_name = config.folder_name
        params = config.params
        self.parameters = params.parameters
        #self.forcing_data = config.forcing_data_file.forcing_data
        self.forcing_data = config.forcing_data

        #self.forcing_data = pd.read_csv(config.forcing_data_string,index_col=0)#,header=0, delim_whitespace=True)
        self.length_of_simulation = config.length_of_simulation

        self.precip_data = np.array(self.forcing_data['P(mm/h)'])[0:self.length_of_simulation] #these could probably happen in init too
        self.PET_data = np.array(self.forcing_data['PET(mm/h)'])[0:self.length_of_simulation]



        self.h_p_init = 0
        self.h_p_vec=[self.h_p_init]
        self.h_p = self.h_p_init
        #self.load_ICs = config.load_ICs

        self.initial_psi = config.initial_psi

        self.parameters = params.parameters

        self.theta_r_vec   = params.theta_r_vec
        self.theta_s_vec   = params.theta_s_vec
        self.K_s_vec       = params.K_s_vec
        self.alpha_vec     = params.alpha_vec
        self.n_vec         = params.n_vec
        self.m_vec         = params.m_vec
        self.max_depth_vec = params.max_depth_vec

        self.h_p_max = params.h_p_max

        # theta_r_vec=self.theta_r_vec #yeah for example don't need these, not super clean
        # theta_s_vec   = params.theta_s_vec
        # K_s_vec       = params.K_s_vec
        # alpha_vec     = params.alpha_vec
        # n_vec         = params.n_vec
        # m_vec         = params.m_vec
        # max_depth_vec = params.max_depth_vec

        self.verbose = config.verbose

        self.parameters = params.parameters
        parameters = self.parameters

        wilting_point_layer_0 = self.theta_of_psi(params.wilting_point_head, params.theta_r_vec[0], params.theta_s_vec[0], params.alpha_vec[0], params.n_vec[0], params.m_vec[0])
        # wilting_point_layer_1 = self.theta_of_psi(params.wilting_point_head, params.theta_r_vec[1], params.theta_s_vec[1], params.alpha_vec[1], params.n_vec[1], params.m_vec[1])
        # wilting_point_layer_2 = self.theta_of_psi(params.wilting_point_head, params.theta_r_vec[2], params.theta_s_vec[2], params.alpha_vec[2], params.n_vec[2], params.m_vec[2])

        #self.wilting_point_vec  = [wilting_point_layer_0,   wilting_point_layer_1,   wilting_point_layer_2]
        self.wilting_point_vec  = [wilting_point_layer_0]



        ###if (self.load_ICs==0):
        theta_of_psi = self.theta_of_psi
        #initial_psi = 20000#100#470
        ###let's say this is the authoritative test for dZdt erroneous boost when crossing layers (this error has been resolved)
        #theta_0 = theta_of_psi(initial_psi, theta_r_vec[0], theta_s_vec[0], alpha_vec[0], n_vec[0], m_vec[0])
        initial_theta_vec = []
        for k in range(0,len(self.parameters)):
            temp_theta_for_init = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[k], params.theta_s_vec[k], params.alpha_vec[k], params.n_vec[k], params.m_vec[k])
            initial_theta_vec.append(temp_theta_for_init)

        # theta_1 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[0], params.theta_s_vec[0], params.alpha_vec[0], params.n_vec[0], params.m_vec[0])#1.6
        # theta_3 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[1], params.theta_s_vec[1], params.alpha_vec[1], params.n_vec[1], params.m_vec[1])#1.9
        # theta_4 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[2], params.theta_s_vec[2], params.alpha_vec[2], params.n_vec[2], params.m_vec[2])#1.9

        #current_states = deque([  [50,theta_0,0], [max_depth_vec[0],theta_1,0], [25,theta_2,1], [max_depth_vec[1],theta_3,1], [max_depth_vec[2],theta_4,2] ])
        self.current_states = deque()
        for b in range(0,len(self.parameters)):
            self.current_states.insert(b,[self.max_depth_vec[b],initial_theta_vec[b],b,b,0])
            #self.current_states = deque([  [params.max_depth_vec[0],theta_1,0,0,0], [params.max_depth_vec[1],theta_3,1,1,0], [params.max_depth_vec[2],theta_4,2,2,0] ])

        #current_states = deque([  [max_depth_vec[0],theta_1,0], [25,theta_2,1], [max_depth_vec[1],theta_3,1], [max_depth_vec[2],theta_5,2] ])#this is a linked list that contains all current state vars, initially set up to contain initial conditions.
        #initial_num_wetting_fronts = 5
        self.initial_num_wetting_fronts = len(self.parameters)
        if self.verbose:
            print(' ')
            print(' ')
            print('Model initialized with the following linked list:')
            print(self.current_states)
        ###else:
        ###    print('loading ICs no longer supported')

        self.states_array = np.array(self.current_states)
        self.wetting_front_numbers = np.array(range(0,self.initial_num_wetting_fronts))

        self.wetting_front_depths = np.array(self.current_states)[:,0]
        self.wetting_front_moistures = np.array(self.current_states)[:,1]

        self.fluxes = [0] #this will record, for each time step, the net fluxes going into the soil
        self.h_p_fluxes = [0]
        self.runoff_vec = [0] #this vector will contain the runoff flux for each time step
        #runoff_vec=np.zeros(length_of_simulation)

        #f_p_vec = [0] ### removed 7 dec

        #f_p_vec = [] #this vector will contain f_p for each time step; if there is no rain, f_p is not calculated and is set to 0, for ease of viewing on a graph (although in reality f_p can be come quite large in periods of no rain). This is not a primary model output and is mostly useful to calcualte actual infiltration, which is of primary interest. Finally, this vector begins with a 0 rather than being empty because at present, the code is not set up to record f_p for the ifrst time step (also implying that the first time step cannot contain precipitation)

        self.actual_infiltration_vec = [0] #this will evventually contain the actual infiltration flux into the soil for each time step

        self.f_p=0

        self.actual_ET_vec = [0]

        self.mass_vec = [0] #this is a vector containing the mass water in the soil profile at each time step; it is used to correct a mass balance error that otherwise would occur if, for a single time step, calculated infiltration would cause a resulting theta value above theta_s for a given layer

        self.free_drainage_flux_vec = [0] #this is the free drainage flux per time step in mm/h

        self.mass_bal_vec = [0]

        self.first_time_step=True
        self.pickle_step = 0


        #loop that runs derivs (dZdt) and then corrects theta via mass balance / psi setting

        self.previous_states = self.current_states

        self.end_time = config.time_step*(config.length_of_simulation)

        self.time_step = config.time_step


    def update(self):
        #updates the model by 1 time step
        self.run_model(self.current_time_step+1)

    def update_until(self,stop_time_step):
        #updates the model until the stop time step; if stop_time_step<self.current_time, then the model won't update.
        self.run_model(stop_time_step)

    def finalize(self):
        print('finalizing model results')
        self.runtime_1 = datetime.now() #this records the end time for computations
        np.save('states_array',self.states_array)

        if (self.first_time_step==False): #ensures that model has to be run for at least 1 time step before it can be saved to disk. It might be that this is not compatible with BMI, and ultimately I mgiht have to change it. I made this choice because specifying fluxes at the start of the model seemserroneous..
            ###this code saves results to the disk; it is more efficient to do this in batches rather than all at once, especially when the run is long.
            #print('len')
            #print(len(forcing_data[(i-100):i]))



            # for i in (range(0,len(states_array))):
            #     mass_bal_vec = np.append(mass_bal_vec, (calc_mass_bal((states_array[states_array[:,4]==i]))))
            # mass_balance_error_vec = mass_bal_vec[0:length_of_simulation]-mass_bal_vec[0]+fluxes_cumulative[0:length_of_simulation]+h_p_vec[0:length_of_simulation]

            #forcing_data_save = self.forcing_data[(i-mod_number):i]
            #forcing_data_save = self.forcing_data[0:self.length_of_simulation]
            forcing_data_save = self.forcing_data[0:self.current_time_step]
            try:
                forcing_data_save['runoff[mm/h]'] = self.runoff_vec
                forcing_data_save['actual_infil[mm/h]'] = self.actual_infiltration_vec
                #h_p_vec = h_p_vec[:len(precip_data)]
                forcing_data_save['ponded_head[mm]'] = self.h_p_vec[0:(len(self.h_p_vec)-1)]
                forcing_data_save['bottom_flux[mm/h]'] = self.free_drainage_flux_vec
                #forcing_data_save['bottom_demand[mm]'] = free_drainage_vec
                forcing_data_save['water_in_soil[mm]'] = self.mass_vec
                forcing_data_save['mass_bal_error(mm)'] = self.mass_bal_vec#mass_balance_error_vec
                forcing_data_save['actual_ET_per_step(mm)'] = self.actual_ET_vec
            except: #not necessary anymore
                forcing_data_save['runoff[mm/h]'] = self.runoff_vec[0:(len(self.runoff_vec)-1)]
                forcing_data_save['actual_infil[mm/h]'] = self.actual_infiltration_vec[0:(len(self.runoff_vec)-1)]
                #h_p_vec = h_p_vec[:len(precip_data)]
                forcing_data_save['ponded_head[mm]'] = self.h_p_vec[0:(len(self.runoff_vec)-1)]
                forcing_data_save['bottom_flux[mm/h]'] = self.free_drainage_flux_vec[0:(len(self.runoff_vec)-1)]
                #forcing_data_save['bottom_demand[mm]'] = free_drainage_vec
                forcing_data_save['water_in_soil[mm]'] = self.mass_vec[0:(len(self.runoff_vec)-1)]
                forcing_data_save['mass_bal_error(mm)'] = self.mass_bal_vec[0:(len(self.runoff_vec)-1)]#mass_balance_error_vec
                forcing_data_save['actual_ET_per_step(mm)'] = self.actual_ET_vec[0:(len(self.runoff_vec)-1)]


            #file_str = "pickle_jar/LGAR_output"+str(self.pickle_step)+".pkl"
            #file_str = self.folder_name + "/LGAR_output"+str(self.pickle_step)+".pkl"
            #file_str = self.folder_name + "/LGAR_output.pkl"
            #file_str = 'pickle_jar' + "/LGAR_output.pkl"

            file_str = self.config.output_file_name_fluxes
            forcing_data_save.to_pickle(file_str)

                # np.save('states_array',states_array)
                #
                # states_array = np.array(current_states)
                # total_num_wetting_fronts = len(current_states)
                # wetting_front_numbers = np.array(range(0,total_num_wetting_fronts))
                #
                # states_array = np.column_stack((np.array(current_states), wetting_front_numbers))#adds wetting front number; code is set up in such a way that each wetting front gets an integer number, starting at 0, from superficial to deep, just to represent the total number of wetting fronts
                #
                # states_array = np.column_stack((states_array, np.zeros(total_num_wetting_fronts)))#adds time step of 0; when current_states is updated later, the time step will be added instead of 0

                # self.pickle_step=self.pickle_step+1
                #
                #
                #
                # #pickle.dump(forcing_data_save,open( "LGAR_output.pkl", "ab" ))
                #
                #
                #
                # self.fluxes=[]
                # #fluxes_cumulative=[]
                # #f_p_vec=[]
                # self.runoff_vec=[]
                # self.actual_infiltration_vec=[]
                # self.h_p_vec=[]
                # self.free_drainage_flux_vec=[]
                # #free_drainage_vec=[]
                # self.mass_vec=[]
                # self.actual_ET_vec=[]
                # self.mass_bal_vec=[]

        print("run success! computation runtime:")
        print(self.runtime_1-self.runtime_0)


    def get_attribute(self, att_name):
        try:
            return self._att_map[ att_name.lower() ]
        except:
            print(' ERROR: Could not find attribute: ' + att_name)

    def get_component_name(self):
        """Name of the component."""
        return self.get_attribute( 'model_name' )

    def get_input_item_count(self):
        """ """
        return len(self._input_var_names)

    def get_output_item_count(self):
        """ """
        return len(self._output_var_names)

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_var_grid(self, name):
        # JG Edit
        # all vars have grid 0 but check if its in names list first
        #PTL: not sure if grid will be 0 for soil moisture; I would prefer that, and report soil moisture in terms of depth-theta pairs that represent wetting fronts
        if name in (self._output_var_names + self._input_var_names):
            return self._var_grid_id

    def get_var_type(self, long_var_name):
        """Data type of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        str
            Data type.
        """
        #return self.get_value_ptr(long_var_name)  #.dtype
        return self._var_types[long_var_name]  #.dtype

    # def get_var_units(self, long_var_name):
    #     return self._var_units[long_var_name]

    def get_var_units(self, long_var_name):
        return self._var_name_units_map[ long_var_name ][1] #_var_units_map will have to be set in initialize

    def get_var_itemsize(self, name): #should probably be long name
#        return np.dtype(self.get_var_type(name)).itemsize
        check = self.get_value(name)
        if check == 'this variable does not yet exist because the model has to be run for at least 1 time step first':
            return(check)
        else:
            return(np.array(self.get_value(name)).itemsize)

    def get_value(self, long_var_name):#, dest
        # name_to_use = self._var_name_units_map[var_name][0]
        # return(self.name_to_use)
        #self._long_names_to_model_variables_map.update()
        #print(getattr( self, self._var_name_units_map[longvar_name][0] ))
        try:
            #return(self._long_names_to_model_variables_map[long_var_name])
            #print(getattr( self, self._var_name_units_map[long_var_name][0] ))
            return( np.array( getattr( self, self._var_name_units_map[long_var_name][0] )).astype('float32') )
        except:
            return('Variable name not recognized. Use get_input_var_names or get_output_var_names for more info.')

        # if var_name == 'current_states':
        #     return self.current_states
        # elif var_name == 'actual_infiltration':
        #     try:
        #         return np.array(self.actual_infiltration_vec[-1])
        #         #exec('self.'+dest+' = np.array(self.actual_infiltration_vec[-1]')
        #     except:
        #         return('this variable does not yet exist because the model has to be run for at least 1 time step first')
        # elif var_name == 'wetting_front_moistures':
        #     try:
        #         return np.array(self.wetting_front_moistures)
        #         #exec('self.'+dest+' = np.array(self.actual_infiltration_vec[-1]')
        #     except:
        #         return('this variable does not yet exist because the model has to be run for at least 1 time step first')
        # else:
        #     return('variable not recognized. Here are the recognized variables:')

    def get_var_nbytes(self, long_var_name):
        """Get units of variable.
        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        Returns
        -------
        int
            Size of data array in bytes.
        """
        # JMFrame NOTE: Had to import sys for this function
        return sys.getsizeof(self.get_value(long_var_name))

    def get_var_location(self, name):
        # JG Edit
        # all vars have location node but check if its in names list first
        if name in (self._output_var_names + self._input_var_names):
            return self._var_loc

    def get_current_time(self):
        return(self.current_time)

    def get_start_time(self):
        return(0.0) #by defintion, at least so far with LGAR

    def get_end_time(self):
        return(self.end_time)

    def get_time_units(self):
        return self.get_attribute( 'time_units' )

    def get_time_step(self):
        return self.get_attribute( 'time_step_size' ) #JG: Edit

    def get_value_ptr(self, long_var_name):
        return(self.get_value(long_var_name))

    def get_value_at_indices(self, var_name, dest, indices):
        #print('TESTING HERE')
        #print( np.array(self.get_value(var_name)).flatten().shape[0] )
        if np.array(self.get_value(var_name)).flatten().shape[0] == 1:
            return self.get_value(var_name)
            #print('TESTING HERE')
        else:
            val_array = np.array(self.get_value(var_name))#.flatten()
            #print(val_array)
            #print(indices)
            #print('TESTING HERE')
            return np.array([val_array[i] for i in indices]).astype('float32')



# self.precip_mm_per_h = precip_data[i]
# self.PET_mm_per_h = PET_data[i]
# self.runoff_mm_per_h = runoff/time_step
# self.actual_infiltration_mm_per_h = precip_mass_to_add/time_step
# self.bottom_flux_mm_per_h = free_drainage_demand/time_step
# self.actual_ET_mm_per_h = actual_ET_demand/time_step
# self.ponding_depth_mm = self.h_0
# self.wetting_front_depths #already defined
# self.wetting_front_moistures #already defined

    def set_value(self, var_name, value):
        setattr( self, self._var_name_units_map[var_name][0], value )

        if var_name == 'atmosphere_water__precipitation_volume_flux':
            self.precip_data[int(self.current_time/self.time_step)] = value

        elif var_name == 'land_surface_water__potential_evapotranspiration_volume_flux':
            self.PET_data[int(self.current_time/self.time_step)] = value

        elif var_name == 'land_surface_water__runoff_volume_flux':
            self.runoff_vec[-1] = value*self.time_step

        elif var_name == 'soil_water__infiltration_volume_flux':
            self.actual_infiltration_vec[-1] = value*self.time_step

        elif var_name == 'soil_water_sat-zone_top__recharge_volume_flux':
            self.free_drainage_flux_vec[-1] = value*self.time_step

        elif var_name == 'land_surface_water__evapotranspiration_volume_flux':
            self.actual_ET_vec[-1] = value*self.time_step

        elif var_name == 'land_surface_water__infiltration_ponding_depth':
            self.h_p_vec[-1] = value

        elif var_name == 'soil_water_wetting-front__depth':
            self.wetting_front_depths = value
            try:
                self.current_states = np.array(self.current_states)
                self.current_states[:,0] = value
            except:
                print('wrong format for setting wetting front depths')

        elif var_name == 'soil_water__volume_fraction':
            self.wetting_front_moistures = value
            try:
                self.current_states = np.array(self.current_states)
                self.current_states[:,1] = value
            except:
                print('wrong format for setting wetting front moistures')

        else:
            print('Variable name not recognized. Use get_input_var_names or get_output_var_names for more info.' )



                                # 'atmosphere_water__precipitation_volume_flux':'float64',
                                # 'land_surface_water__potential_evapotranspiration_volume_flux':'float64',
                                # 'land_surface_water__runoff_volume_flux':'float64',
                                # 'soil_water__infiltration_volume_flux':'float64',
                                # 'soil_water_sat-zone_top__recharge_volume_flux':'float64',
                                # 'land_surface_water__evapotranspiration_volume_flux':'float64',
                                # 'land_surface_water__infiltration_ponding_depth':'float64',
                                # 'soil_water_wetting-front__depth':'float64', #this will be a vector containing at least 3 float64 numbers
                                # 'soil_water__volume_fraction':'float64' #this will be a vector containing at least 3 float64 numbers


    def set_value_at_indices(self, name, inds, src):
        if np.array(self.get_value(name)).flatten().shape[0] == 1:
            self.set_value(name, src)
        else:
            val = self.get_value_ptr(name)
            val = np.array(val) #because if you want to do this with states_array, it's originally a deque object which has no flatten attribute (method?)
            val = val.flatten()
            # print('TEST VAL HERE')
            # print(val)
            # print(inds.shape)
            for i in inds:
                # print('index i')
                # print(i)
                # print(val.flatten()[inds[i]])
                # print(src[i])
                # print([inds[i]])
                val[i] = src[i]
                #val.flatten()[inds[i-1]] = src[i-1]
            # print('THING THAT WILL BE SET')
            # print(val)
            self.set_value(name, val)







    def get_grid_edge_count(self, grid):
        raise NotImplementedError("get_grid_edge_count")

    #------------------------------------------------------------
    def get_grid_edge_nodes(self, grid, edge_nodes):
        raise NotImplementedError("get_grid_edge_nodes")

    #------------------------------------------------------------
    def get_grid_face_count(self, grid):
        raise NotImplementedError("get_grid_face_count")

    #------------------------------------------------------------
    def get_grid_face_edges(self, grid, face_edges):
        raise NotImplementedError("get_grid_face_edges")

    #------------------------------------------------------------
    def get_grid_face_nodes(self, grid, face_nodes):
        raise NotImplementedError("get_grid_face_nodes")

    #------------------------------------------------------------
    def get_grid_node_count(self, grid):
        raise NotImplementedError("get_grid_node_count")

    #------------------------------------------------------------
    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        raise NotImplementedError("get_grid_nodes_per_face")

    #------------------------------------------------------------
    def get_grid_origin(self, grid_id, origin):
        raise NotImplementedError("get_grid_origin")

    #------------------------------------------------------------
    def get_grid_rank(self, grid_id):

        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return 1

    #------------------------------------------------------------
    def get_grid_shape(self, grid_id, shape):
        raise NotImplementedError("get_grid_shape")

    #------------------------------------------------------------
    def get_grid_size(self, grid_id):

        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return 1

    #------------------------------------------------------------
    def get_grid_spacing(self, grid_id, spacing):
        raise NotImplementedError("get_grid_spacing")

    #------------------------------------------------------------
    def get_grid_type(self, grid_id=0):

        # JG Edit
        # 0 is the only id we have
        if grid_id == 0:
            return 'scalar'

    #------------------------------------------------------------
    def get_grid_x(self):
        raise NotImplementedError("get_grid_x")

    #------------------------------------------------------------
    def get_grid_y(self):
        raise NotImplementedError("get_grid_y")

    #------------------------------------------------------------
    def get_grid_z(self):
        raise NotImplementedError("get_grid_z")
