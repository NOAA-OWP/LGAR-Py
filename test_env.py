#python test environment for LGAR

import importlib
import numpy as np

#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_Panama.py'
config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_synth_1.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_synth_2.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_synth_3.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_Phillipsburg.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_Bushland.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_files/config_Fort_Assiniboine.py'

# folder_name = 'Bushland_project'
# config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/Bushland_project/config.py'

# folder_name = 'Phillipsburg_project'
# config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/Phillipsburg_project/config.py'

# folder_name = 'Fort_Assiniboine_project'
# config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/Fort_Assiniboine_project/config.py'

#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_Phillipsburg.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_synth_1.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_synth_2.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_synth_3.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_synth_4.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_Bushland.py'



#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_Los_Lunas.py'
#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_Sidney.py'

#####config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_video.py'

#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_video_old.py'

#config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/config_short_synthetic.py'


# folder_name = 'Panama_project'
# config_string = '/Users/peterlafollette/Desktop/alt-modular/Modules/GAR/LGAR_BMI/Panama_project/config.py'


#config = importlib.import_module('config.py')

import LGAR_compute as LGAR_compute

test=LGAR_compute.LGAR(config_string)

# print('set_value test')
# test_val=1
# test.set_value('land_surface_water__infiltration_ponding_depth',test_val)
# print(test.actual_infiltration_vec[-1])
# print(test.actual_infiltration_mm_per_h)
# print(' ')
#
# 1/0

test.run_model(test.length_of_simulation)
#test.finalize()
# print('this text means that the class ran successfully')
# print(' ')

test.finalize()



###########everything below here tests the rest of the BMI functions

# print('get_value test, soil moisture vector')
# print(test.get_value('soil_water__volume_fraction'))
# print(' ')
#
# print('update test')
# print(test.update())
# print(' ')
#
# print('update test')
# print(test.update())
# print(' ')
#
# print('update test')
# print(test.update())
# print(' ')
#
# # print('before update')
# # print(test.current_states)
# #test.update_until(3)
# # test.update()
# # print('after update')
# # print(test.current_states)
#
#
# print('finalization test')
# test.finalize()
# print(' ')
#
# # print('init test')
# # test.initialize(config_string)
# # print(test.current_states)
# # print(' ')
#
# print('CURRENT TEST')
# print('get_value test')
# print(test.get_value('soil_water__infiltration_volume_flux'))
# print(' ')
#
#
# print('get_component_name test')
# print(test.get_component_name())
# print(' ')
#
# print('get_input_item_count test')
# print(test.get_input_item_count())
# print(' ')
#
# print('get_output_item_count test')
# print(test.get_output_item_count())
# print(' ')
#
# print('get_input_var_names test')
# print(test.get_input_var_names())
# print(' ')
#
# print('get_output_var_names test')
# print(test.get_output_var_names())
# print(' ')
#
# print('get_var_grid test')
# print(test.get_var_grid('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_var_type test')
# print(test.get_var_type('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_var_units test')
# print(test.get_var_units('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_var_itemsize test') #this should use the long name
# print(test.get_var_itemsize('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_var_nbytes test')
# print(test.get_var_nbytes('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_var_location test')
# print(test.get_var_location('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_current_time test')
# print(test.get_current_time())
# print(' ')
#
# print('get_start_time test')
# print(test.get_start_time())
# print(' ')
#
# print('get_end_time test')
# print(test.get_end_time())
# print(' ')
#
# print('get_time_units test')
# print(test.get_time_units())
# print(' ')
#
# print('get_time_step test')
# print(test.get_time_step())
# print(' ')
#
# print('get_value test, infiltration')
# print(test.get_value('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_value test, soil moisture vector')
# print(test.get_value('soil_water__volume_fraction'))
# print(' ')
#
# print('get_value_ptr test')
# print(test.get_value_ptr('soil_water__infiltration_volume_flux'))
# print(' ')
#
# print('get_value_at_indices test, infiltration')
# print(test.get_value_at_indices('soil_water__infiltration_volume_flux','foo',1))
# print(' ')
#
# print('get_value_at_indices test, vector of moistures')
# print(test.get_value_at_indices('soil_water__volume_fraction','foo',np.array([1,2])))
# print(' ')
#
# print(test.wetting_front_depths)
#
# print('get_value_at_indices test, current_states (should fail)')
# print(test.get_value_at_indices('current_states','foo',np.array([1,2,3])))
# print(' ')
#
# print('set_value test')
# test_val=1
# test.set_value('soil_water__infiltration_volume_flux',test_val)
# print(test.actual_infiltration_vec[-1])
# print(test.actual_infiltration_mm_per_h)
# print(' ')
#
# print('get_value test, infiltration (checking that set value works)')
# print(test.get_value('soil_water__infiltration_volume_flux'))
# print(' ')
#
#
# print('set_value_at_indices test')
# print(test.set_value_at_indices('soil_water__infiltration_volume_flux',np.array([0]),10))
# print(test.actual_infiltration_vec[-1])
# print(test.actual_infiltration_mm_per_h)
# print(' ')
#
# # print('set_value_at_indices test')
# # print(test.set_value_at_indices('wetting_front_moistures',np.array([0]),np.array([0.4])))
# # print(' ')
# # print('using get_value to test set_value_at_indices')
# # print(test.get_value('wetting_front_moistures'))
# # print(' ')
#
# print('set_value_at_indices test')
# print(test.set_value_at_indices('soil_water__volume_fraction',np.array([0,1,2]),np.array([0.4,0.3,0.2])))
# print(' ')
# print('using get_value to test set_value_at_indices')
# print(test.get_value('soil_water__volume_fraction'))
# print(' ')
#
# print('init test')
# test.initialize(config_string)
# print(test.current_states)
# print(' ')
#


# print('set_value_at_indices test')
# print(test.set_value_at_indices('current_states',np.array([0,1,2]),np.array([-9999,-9999,-9999])))
# print(' ')
# print('using get_value to test set_value_at_indices')
# print(test.get_value('current_states'))
# print(' ')

# print('get_grid_edge_count test')
# print(test.get_grid_edge_count('no_grid'))
# print(' ')





# test.update_until(1500)
#
# test.finalize()










#####trying to make my code gel with BMI unit test

# from pathlib import Path
# cfg_file=Path(config_string)
# print(cfg_file)
# test.initialize(cfg_file)
#
#
# print(' ')
# print(test.get_output_var_names()[0])
# var_name = test.get_output_var_names()[0]
# print(' ')
#
# this_set_value = -99.0
# #var_name = 'atmosphere_water__precipitation_volume_flux'
# test.set_value(var_name, this_set_value)
# print ("  set value: " + str(this_set_value))



# 'atmosphere_water__precipitation_volume_flux':'float32',
# 'land_surface_water__potential_evapotranspiration_volume_flux':'float32',
# 'land_surface_water__runoff_volume_flux':'float32',
# 'soil_water__infiltration_volume_flux':'float32',
# 'soil_water_sat-zone_top__recharge_volume_flux':'float32',
# 'land_surface_water__evapotranspiration_volume_flux':'float32',
# 'land_surface_water__infiltration_ponding_depth':'float32',
# 'soil_water_wetting-front__depth':'float32', #this will be a vector containing at least 3 float64 numbers
# 'soil_water__volume_fraction':'float32' #thi
