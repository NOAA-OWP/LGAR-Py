import numpy as np

#here, layer parameters are defined.

num_layers = 3

#first layer: loam (HYDRUS catalog)
theta_r_layer_0          = 0.0416  #residual volumetric water content.            unitless
theta_s_layer_0          = 0.4189   #total porosity.                               unitless
K_s_layer_0              = 10.2   #saturated hydraulic conductivity.             mm/h
###fine-coarse-fine testing
#K_s_layer_0              = 1.04   #saturated hydraulic conductivity.             mm/h
alpha_layer_0            = 0.002393 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_0                = 1.3527   #unitless
m_layer_0                = 1-1/n_layer_0   #unitless
max_depth_layer_0        = 130 #100   #mm, this is the thickness of the first soil layer
#note that alpha, n, and m are van Genuchten parameters
param_set_0 = [theta_r_layer_0, theta_s_layer_0, K_s_layer_0, alpha_layer_0, n_layer_0, m_layer_0, max_depth_layer_0]


#second layer: clay loam (HYDRUS catalog)
theta_r_layer_1         = 0.0762  #residual volumetric water content.            unitless
theta_s_layer_1         = 0.4794   #total porosity.                               unitless
K_s_layer_1             = 2.6    #saturated hydraulic conductivity.             mm/h
alpha_layer_1           = 0.001648 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_1               = 1.3086   #unitless
m_layer_1               = 1-1/n_layer_1   #unitless
max_depth_layer_1       = 710    #mm
#note that alpha, n, and m are van Genuchten parameters
param_set_1 = [theta_r_layer_1, theta_s_layer_1, K_s_layer_1, alpha_layer_1, n_layer_1, m_layer_1, max_depth_layer_1]



#second layer: silty clay loam (finer than layer above) (HYDRUS catalog)
theta_r_layer_2         = 0.0574  #residual volumetric water content.            unitless
theta_s_layer_2         = 0.3986   #total porosity.                               unitless
K_s_layer_2             = 2.6    #saturated hydraulic conductivity.             mm/h
alpha_layer_2           = 0.000458 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_2               = 1.4243   #unitless
m_layer_2               = 1-1/n_layer_2   #unitless
max_depth_layer_2       = 1700 #10000   #mm
#note that alpha, n, and m are van Genuchten parameters
param_set_2 = [theta_r_layer_2, theta_s_layer_2, K_s_layer_2, alpha_layer_2, n_layer_2, m_layer_2, max_depth_layer_2]


#also making vectors of parameters based on parameter name
theta_r_vec   = [theta_r_layer_0,   theta_r_layer_1,   theta_r_layer_2]
theta_s_vec   = [theta_s_layer_0,   theta_s_layer_1,   theta_s_layer_2]
K_s_vec       = [K_s_layer_0,       K_s_layer_1,       K_s_layer_2]
alpha_vec     = [alpha_layer_0,     alpha_layer_1,     alpha_layer_2]
n_vec         = [n_layer_0,         n_layer_1,         n_layer_2]
m_vec         = [m_layer_0,         m_layer_1,         m_layer_2]
max_depth_vec = [max_depth_layer_0, max_depth_layer_1, max_depth_layer_2]

#here, the parameter h_p_max, or the maximum ponding depth (in mm), is defined. It is defined separately because it is the only parameter that is independent of individual layers.
h_p_max = 0 #20

wilting_point_head = 154950 #this is the pressure head in mm (but negative) that represents the pressure of the wilting point. A value of 154950 mm physically indicates, in LGAR, a wilting point of -15 atmospheres. This is used to calculate the volumetric water content which represents the wilting point.
relative_moisture_at_which_PET_equals_AET = 0.75 #0.75 #technically as currently implemented this is just the factor by which theta_s is multiplied





###
###for testing what happens when 2 adjacent layers are the same
# theta_r_vec   = [theta_r_layer_0,   theta_r_layer_0,   theta_r_layer_2]
# theta_s_vec   = [theta_s_layer_0,   theta_s_layer_0,   theta_s_layer_2]
# K_s_vec       = [K_s_layer_0,       K_s_layer_0,       K_s_layer_2]
# alpha_vec     = [alpha_layer_0,     alpha_layer_0,     alpha_layer_2]
# n_vec         = [n_layer_0,         n_layer_0,         n_layer_2]
# m_vec         = [m_layer_0,         m_layer_0,         m_layer_2]
# max_depth_vec = [max_depth_layer_0, max_depth_layer_0, max_depth_layer_2]
# param_set_1 = param_set_0
###

#also putting parameter sets into array because it is necessary for loop that considers all wetting fronts in all layers
parameters = np.vstack([param_set_0,param_set_1,param_set_2])

#this will be useful later; just defining what the number of the maximum possible layer is
max_layer = len(parameters)-1


#finally, these few lines establish what the exact boundary depths between layers are, which will be useful later for propagating wetting fronts between layers.
boundary_depths = []
temp_bdy_depth = 0
for s in range(0,len(parameters)):
    temp_bdy_depth = temp_bdy_depth + parameters[s][-1]
    boundary_depths.append(temp_bdy_depth)
