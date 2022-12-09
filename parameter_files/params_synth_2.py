import numpy as np

###In the first block, layer parameters are defined. The second block does not need to be adjusted by the user.



#######

#here, the parameter h_p_max, or the maximum ponding depth (in mm), is defined.
h_p_max = 0
###the following two parameters relate to the correction of PET to AET.
wilting_point_head = 154950 #this is the pressure head in mm (but negative) that represents the pressure of the wilting point. A value of 154950 mm physically indicates, in LGAR, a wilting point of -15 atmospheres. This is used to calculate the volumetric water content which represents the wilting point in the top layer.
relative_moisture_at_which_PET_equals_AET = 0.75 #relative moisture at which AET is approximately the same as PET
num_layers = 3 #total number of soil layers that are in the vadose zone

##first layer: silt (HYDRUS catalog)
theta_r_layer_0          = 0.034  #residual volumetric water content.            unitless
theta_s_layer_0          = 0.46#43   #total porosity.                               unitless
K_s_layer_0              = 2.5  #saturated hydraulic conductivity.             mm/h
###fine-coarse-fine testing
#K_s_layer_0              = 1.04   #saturated hydraulic conductivity.             mm/h
alpha_layer_0            = 0.0016 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_0                = 1.37   #unitless
m_layer_0                = 1-1/n_layer_0   #unitless
max_depth_layer_0        = 100    #mm, this is the thickness of the first soil layer
#note that alpha, n, and m are van Genuchten parameters
param_set_0 = [theta_r_layer_0, theta_s_layer_0, K_s_layer_0, alpha_layer_0, n_layer_0, m_layer_0, max_depth_layer_0]

#second layer: silty clay loam (HYDRUS catalog)
theta_r_layer_1         = 0.089  #residual volumetric water content.            unitless
theta_s_layer_1         = 0.43   #total porosity.                               unitless
K_s_layer_1             = 0.7    #saturated hydraulic conductivity.             mm/h
alpha_layer_1           = 0.0010 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_1               = 1.23   #unitless
m_layer_1               = 1-1/n_layer_1   #unitless
max_depth_layer_1       = 300    #mm
#note that alpha, n, and m are van Genuchten parameters
param_set_1 = [theta_r_layer_1, theta_s_layer_1, K_s_layer_1, alpha_layer_1, n_layer_1, m_layer_1, max_depth_layer_1]

#third layer: (HYDRUS catalog)
theta_r_layer_2         = 0.068  #residual volumetric water content.            unitless
theta_s_layer_2         = 0.38 #43  #total porosity.                               unitless
K_s_layer_2             = 2      #saturated hydraulic conductivity.             mm/h
alpha_layer_2           = 0.0008 #keep in mind the unit here is mm^-1, often alpha is reported in cm^-1
n_layer_2               = 1.09   #unitless
m_layer_2               = 1-1/n_layer_2   #unitless
max_depth_layer_2       = 300    #mm
#note that alpha, n, and m are van Genuchten parameters
param_set_2 = [theta_r_layer_2, theta_s_layer_2, K_s_layer_2, alpha_layer_2, n_layer_2, m_layer_2, max_depth_layer_2]

#######









#######

#also putting parameter sets into array because it is necessary for loop that considers all wetting fronts in all layers
parameters = np.vstack([param_set_0,param_set_1,param_set_2])

if (len(parameters)>num_layers):
    while (len(parameters)>num_layers):
        parameters = np.delete(parameters, len(parameters)-1, 0)

theta_r_vec = []
theta_s_vec = []
K_s_vec = []
alpha_vec = []
n_vec = []
m_vec = []
max_depth_vec = []

for param_set_num in range(0,num_layers):
    theta_r_vec.append(parameters[param_set_num][0])
    theta_s_vec.append(parameters[param_set_num][1])
    K_s_vec.append(parameters[param_set_num][2])
    alpha_vec.append(parameters[param_set_num][3])
    n_vec.append(parameters[param_set_num][4])
    m_vec.append(parameters[param_set_num][5])
    max_depth_vec.append(parameters[param_set_num][6])


#this will be useful later; just defining what the number of the maximum possible layer is
max_layer = len(parameters)-1

#finally, these few lines establish what the exact boundary depths between layers are, which will be useful later for propagating wetting fronts between layers.
boundary_depths = []
temp_bdy_depth = 0
for s in range(0,len(parameters)):
    temp_bdy_depth = temp_bdy_depth + parameters[s][-1]
    boundary_depths.append(temp_bdy_depth)

#######
