

### this fairly large .py file does the computations for LGAR. the LGAR model can be run by typing "python LGAR_compute.py" in terminal, once in the correct folder (containing all the files from the LGAR_standalone repo).
### Currently, this file will run for a number of time steps specified on the second to last line of forcing_data.py. Meaning, it currently has a sort of "run until" capability, and actually if the specified runtime of the simulation was just 1 time step, then I suppose it also has a sort of "advance by 1 time step" functionality too.

from collections import deque
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import copy
import cProfile
from pathlib import Path

import BMI_base_class as BMI


# config_string = "config"
# command_string = "import " + config_string + " as config"
# exec(command_string)
#
# import config as config
# params = config.params
# parameters = params.parameters
# forcing_data = config.forcing_data_file.forcing_data
# length_of_simulation = config.length_of_simulation
#
# precip_data = np.array(forcing_data['P(mm/h)'])[0:length_of_simulation] #these could probably happen in init too
# PET_data = np.array(forcing_data['PET(mm/h)'])[0:length_of_simulation]







#intermediate (helper) fxns, which help us convert between psi and theta and capital theta and calculate capillary drive
#these functions can be used for any layer, as long as the correct parameters are passed for that layer.
#for example, if theta_now is equal to the theta value of a wetting front in the first layer, then
#psi for that wetting front could be calculated with
#psi_of_theta(theta_now,theta_s_layer_0,theta_r_layer_0,n_layer_0,m_layer_0,alpha_layer_0)
#and note that these functions take variables as arguments first, then parameters

#psi values in LGAR are either 0 or positive. However, LGAR (and seemingly most or all Green-Ampt based models) are formualted using the absolute value of psi. So for example, when a pis value of 100 mm is used in LGAR, this physically indicates unsaturated soil with a psi value of -100 mm. A psi value of 0 still represents saturation.

def capital_theta(theta, theta_r, theta_s): #relative water content, such that this yields 1 [unitless] if theta=theta_s (saturation) and 0 if theta=theta_r (theta_r is the residual water content)
    return( (theta - theta_r) / (theta_s - theta_r) )

def K(capital_theta_temp,K_s,m): #hydraulic conductivity (mm/h) as a function of capital theta, taken from some of Fred's code (a file called one_block)
    return(  K_s*(capital_theta_temp)**0.5 * (1-(1-(capital_theta_temp)**(1/m))**m)**2  )

def psi_of_theta(theta,theta_s,theta_r,n,m,alpha): #gives the hydraulic head [mm], given the volumetric water content and some other parameters
    return(  ((((theta_s-theta_r)/(theta-theta_r))**(1/m)-1)**(1/n))/alpha  )

def capital_theta_of_psi(psi,alpha,n,m): #relative water content [unitless], based on hydraulic head
    return (  1/((1+(alpha*psi)**n)**m)  )

#capillary drive function (final units = mm), adapted from equation 8 in 1997 GAR paper "Green and Ampt Infiltration with Redistribution", Ogden and Saghafian
#this equation integrates the hydraulic conductivity over psi (hydraulic head) to determine how strong the capillary suction is
#it uses the trapezoidal rule for integration with a default number of trapezoids of 20.
#using a smaller number will result in a less precise value for G but will improve computation time. So far testing has suggested that you shouldn't go below 20.
#finally as of ~February 2022, Scott Peckham has been working on a closed form for G. If successful, this should take ~30% off the runtime of LGAR.


def G_closed(theta_upper_lim,theta_lower_lim,psi_b,lambdaa,theta_r,theta_s):
    H_c = psi_b*((2+3*lambdaa)/(1+3*lambdaa))
    capital_theta_current = capital_theta(theta_upper_lim,theta_r,theta_s)
    capital_theta_below = capital_theta(theta_lower_lim,theta_r,theta_s)
    try:
        G_result = H_c*(capital_theta_current**(3+1/lambdaa)-capital_theta_below**(3+1/lambdaa))/(1-capital_theta_below**(3+1/lambdaa))
    except:
        G_result = H_c
    return(G_result)



def G(psi_upper_lim, psi_lower_lim, alpha, n, m, K_s, num_trapezoids=120):#120 seems to be a good value for num_trapezoids
    if (psi_upper_lim==psi_lower_lim):
        return(0)
    else:
        ###multiple experiments were done here to try to optimize the speed of integration. The uncommented code was fastest. Uncommenting the following 2 lines (and commenting out the currently active psis_to_integrate line) will also work but is slower. Keeping it in as comments however.
        # psis_to_integrate=np.arange(psi_lower_lim, psi_upper_lim, ((psi_upper_lim-psi_lower_lim)/num_trapezoids)).tolist()
        # psis_to_integrate.append(psi_upper_lim)
        psis_to_integrate=np.arange(psi_lower_lim, psi_upper_lim, ((psi_upper_lim-psi_lower_lim)/num_trapezoids+1e-10))
        #print(len(psis_to_integrate))
        #print(psis_to_integrate)
        #psis_to_integrate=np.arange(psi_lower_lim, psi_upper_lim, ((psi_upper_lim-psi_lower_lim)/num_trapezoids)).tolist()
        #psis_to_integrate.append(psi_upper_lim)
        # print(len(psis_to_integrate))
        # print(psis_to_integrate)
        # print('psi limits')
        # print(psi_lower_lim)
        # print(psi_upper_lim)

        G_result=0
        for i in range(0,num_trapezoids):
            G_result=G_result+(K(capital_theta_of_psi(psis_to_integrate[i+1],alpha,n,m),K_s,m)+K(capital_theta_of_psi(psis_to_integrate[i],alpha,n,m),K_s,m))/2*abs((psis_to_integrate[i+1]-psis_to_integrate[i]))
        G_result=abs(G_result/K_s)
        if (psi_upper_lim>psi_lower_lim):
            G_result = abs(G_result) #technically, if the layer below is wetter than the layer above, then you will have upward capillary suction.
            #however, in the LGAR concept, we only have increasing psi with depth. Allowing for drier wetting fronts within the same layer is a possible advancement of LGAR (and this is probably something we want in LGARTO)
        return(G_result)

def theta_of_psi(psi, theta_r, theta_s, alpha, n, m): #gives volumetric water content (unitless) based on hydraulic head
    # print('hey look here here here is psi')
    # print(psi)
    return(  theta_r + (theta_s-theta_r)/(1+(alpha*psi)**n)**m  )

def calc_abs_depth(deque_row): #later, depths of wetting fronts will be defined as the depth of a wetting front in its layer, rather than its absolute depth (for example, a wetting front that is 25 mm deep into layer 1, when layer 0 is 100 mm thick, has an absolute depth of 125 mm). It is however sometimes necessary to know the absolute depth of a wetting front rather than just its depth in a layer
    current_layer = deque_row[2]
    abs_depth = deque_row[0]
    for q in range(0,current_layer):
        abs_depth = abs_depth + parameters[q][-1]
    return(abs_depth)

def using_list_comp(sample, indices):
    return tuple([t[i] for t in sample] for i in indices)

#technically, this function just calculates the current mass of water in the soil; it does not actually calculate the mass balance error, which is done later.
def calc_mass_bal(current_states):
    mass = 0
    for l in (range(0,len(current_states))):
        if (l<(len(current_states)-1)):
            current_layer = current_states[l][2]
            next_layer = current_states[l+1][2]
            if (next_layer==current_layer):
                temp_mass = current_states[l][0]*(current_states[l][1]-current_states[l+1][1])
                mass = mass + temp_mass
            else:
                temp_mass = (current_states[l][0])*current_states[l][1]
                mass = mass + temp_mass
        if l==(len(current_states)-1):
            temp_mass = (current_states[l][0])*current_states[l][1]
            mass = mass + temp_mass

    return(mass)






#the function derivs is defined. It takes as an input current_states, which is a linked list (deque object) that contains depth and theta info for all wetting fronts (and more recently edited to contain "wetting front number", which is the number of a wetting front from superficial to deep, and the time step).
#it returns derivatives for the depth of every wetting front, dZdt, and a value of 0 for dthetadt. Theta is later updated by mass balance, as each wetting front conserves its mass.






class LGAR(BMI.base_BMI):



    def __init__(self,config_string):#can add config string as argument
        super(LGAR, self).__init__() #this makes it so I get the stuff from __init__ in the base class as well rather than just overwriting the __init__ function
        self.initialize(config_string)



##### testing if __init__ can just call initialize
        # config_string = "config"
        # command_string = "import " + config_string + " as config"
        # exec(command_string)
        # print("TEST OF IMPORTIING IN INIT")
        # print(config.params)
        #
        # params = config.params
        #
        # self.h_p_init = 0
        # self.h_p_vec=[]
        # self.h_p = self.h_p_init
        # self.load_ICs = config.load_ICs
        #
        # self.initial_psi = config.initial_psi
        #
        # self.parameters = params.parameters
        #
        # self.theta_r_vec   = params.theta_r_vec
        # self.theta_s_vec   = params.theta_s_vec
        # self.K_s_vec       = params.K_s_vec
        # self.alpha_vec     = params.alpha_vec
        # self.n_vec         = params.n_vec
        # self.m_vec         = params.m_vec
        # self.max_depth_vec = params.max_depth_vec
        #
        # # theta_r_vec=self.theta_r_vec #yeah for example don't need these, not super clean
        # # theta_s_vec   = params.theta_s_vec
        # # K_s_vec       = params.K_s_vec
        # # alpha_vec     = params.alpha_vec
        # # n_vec         = params.n_vec
        # # m_vec         = params.m_vec
        # # max_depth_vec = params.max_depth_vec
        #
        # self.verbose = config.verbose
        #
        # self.parameters = params.parameters
        # parameters = self.parameters
        #
        # if (self.load_ICs==0):
        #
        #     #initial_psi = 20000#100#470
        #     ###let's say this is the authoritative test for dZdt erroneous boost when crossing layers (this error has been resolved)
        #     #theta_0 = theta_of_psi(initial_psi, theta_r_vec[0], theta_s_vec[0], alpha_vec[0], n_vec[0], m_vec[0])
        #     theta_1 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[0], params.theta_s_vec[0], params.alpha_vec[0], params.n_vec[0], params.m_vec[0])#1.6
        #     #theta_2 = theta_of_psi(initial_psi*1.0, theta_r_vec[1], theta_s_vec[1], alpha_vec[1], n_vec[1], m_vec[1])#1.6
        #     theta_3 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[1], params.theta_s_vec[1], params.alpha_vec[1], params.n_vec[1], params.m_vec[1])#1.9
        #     #theta_3 = theta_of_psi(initial_psi*1.9, theta_r_vec[1], theta_s_vec[1], alpha_vec[1], n_vec[1], m_vec[1])
        #     #theta_4 = theta_of_psi(initial_psi*1.9, theta_r_vec[2], theta_s_vec[2], alpha_vec[2], n_vec[2], m_vec[2])
        #     theta_4 = theta_of_psi(config.initial_psi*1.0, params.theta_r_vec[2], params.theta_s_vec[2], params.alpha_vec[2], params.n_vec[2], params.m_vec[2])#1.9
        #
        #     #current_states = deque([  [50,theta_0,0], [max_depth_vec[0],theta_1,0], [25,theta_2,1], [max_depth_vec[1],theta_3,1], [max_depth_vec[2],theta_4,2] ])
        #     self.current_states = deque([  [params.max_depth_vec[0],theta_1,0,0,0], [params.max_depth_vec[1],theta_3,1,1,0], [params.max_depth_vec[2],theta_4,2,2,0] ])
        #     #current_states = deque([  [max_depth_vec[0],theta_1,0], [25,theta_2,1], [max_depth_vec[1],theta_3,1], [max_depth_vec[2],theta_5,2] ])#this is a linked list that contains all current state vars, initially set up to contain initial conditions.
        #     #initial_num_wetting_fronts = 5
        #     self.initial_num_wetting_fronts = 3
        # else:
        #     print('loading ICs no longer supported')
        #
        # self.states_array = np.array(self.current_states)
        # self.wetting_front_numbers = np.array(range(0,self.initial_num_wetting_fronts))
        #
        #
        #
        # self.fluxes = [] #this will record, for each time step, the net fluxes going into the soil
        #
        # self.runoff_vec = [] #this vector will contain the runoff flux for each time step
        # #runoff_vec=np.zeros(length_of_simulation)
        #
        # #f_p_vec = [0] ### removed 7 dec
        #
        # #f_p_vec = [] #this vector will contain f_p for each time step; if there is no rain, f_p is not calculated and is set to 0, for ease of viewing on a graph (although in reality f_p can be come quite large in periods of no rain). This is not a primary model output and is mostly useful to calcualte actual infiltration, which is of primary interest. Finally, this vector begins with a 0 rather than being empty because at present, the code is not set up to record f_p for the ifrst time step (also implying that the first time step cannot contain precipitation)
        #
        # self.actual_infiltration_vec = [] #this will evventually contain the actual infiltration flux into the soil for each time step
        #
        # self.f_p=0
        #
        # self.actual_ET_vec = []
        #
        # self.mass_vec = [] #this is a vector containing the mass water in the soil profile at each time step; it is used to correct a mass balance error that otherwise would occur if, for a single time step, calculated infiltration would cause a resulting theta value above theta_s for a given layer
        #
        # self.free_drainage_flux_vec = [] #this is the free drainage flux per time step in mm/h
        #
        # self.mass_bal_vec = []
        #
        # self.first_time_step=True
        # self.pickle_step = 0
        #
        #
        # #loop that runs derivs (dZdt) and then corrects theta via mass balance / psi setting
        #
        # self.previous_states = self.current_states
#####






    #here initial conditions are set. In most cases, setting realistic initial conditions is typically established via a spin up period, although initial conditions can be specified here.

    #This model runs via updating an object called current_states, which is a linked list, for every time step. After current_states is updated, it is appended to states_array.
    #note the syntax for making linked lists with deque, e.g.:
    #current_states = deque([ [30,theta_0,0], [100,theta_1,0], [25,theta_2,1], [100,theta_3,1], [25,theta_4,2], [100,theta_5,2] ])#this is a linked list that contains all current state vars, initially set up to contain initial conditions.
    #this will initalize LGAR with 6 wetting fronts, each indicated by a triplet between brackets, and arranged from superficial to deep.
    #the triplets come in the form [Z,theta,layer_number], where Z is the wetting front depth in its layer, theta is the wetting front moisture (volumetric water content), and layer_number is the layer the wetting front is in.
    #also, LGAR requires that wetting fronts that are adjacent but do not share the same layer have the same psi value.

    #I currently have a lot of different ICs to test different things, although the last one that is not commented out in this cell is the one that is currently used.

    #first, initial conditions in h_p are defined

    #3 feb takeout

    #h_p_vec = [h_p_init]

    #initial conditions



        #from ICs import *



    def calc_abs_depth(self,deque_row): #later, depths of wetting fronts will be defined as the depth of a wetting front in its layer, rather than its absolute depth (for example, a wetting front that is 25 mm deep into layer 1, when layer 0 is 100 mm thick, has an absolute depth of 125 mm). It is however sometimes necessary to know the absolute depth of a wetting front rather than just its depth in a layer
        current_layer = int(deque_row[2])
        abs_depth = deque_row[0]
        for q in range(0,current_layer):
            abs_depth = abs_depth + self.parameters[q][-1]
        return(abs_depth)

    def theta_of_psi(psi, theta_r, theta_s, alpha, n, m): #gives volumetric water content (unitless) based on hydraulic head
        return(  theta_r + (theta_s-theta_r)/(1+(alpha*psi)**n)**m  )


    def derivs (self,current_states,h_p,wf_number_for_h_p):

        config = self.config
        params = config.params
        parameters = params.parameters

        theta_r_vec   = params.theta_r_vec
        theta_s_vec   = params.theta_s_vec
        K_s_vec       = params.K_s_vec
        alpha_vec     = params.alpha_vec
        n_vec         = params.n_vec
        m_vec         = params.m_vec
        max_depth_vec = params.max_depth_vec
        ###again, testing for speed here was done. It turns out that using a list (like f=[]) is the same speed as establishing fixed length array of zeros.
        #f=np.zeros(len(current_states)) #first, f, or the complete array of derivatives for each wetting front, is defined. Again, dZdt is calculated for each wetting front, and dthetadt is set to 0 for each wetting front (theta will later be updated via mass balance).
        f=[]
        #print(f)
        self.number_of_wetting_fronts = len(self.current_states) #because current states is a linked list containing all wetting fronts, its length is the number of wetting fronts currently existing
        maximum_wetting_front_number = self.number_of_wetting_fronts-1 #the wetting fronts are numbered from 0 to n, where wetting front 0 is the most superficial and wetting front n is the deepest. Because this index starts from 0, the number of the deepest wetting front is the length of current-states -1.
        h_p_store = h_p
        layers_vec = using_list_comp(self.current_states,[2])[0]
        for l in range(0,self.number_of_wetting_fronts): #in the following, l will be the wetting front number; derivatives for Z for each wetting front are calculated here.
            Z, theta, layer_number, number_in_soil = self.current_states[l][0:4]

            try:
                Z_below_deriv, theta_below_deriv, layer_number_below_deriv, number_in_soil_below_deriv = self.current_states[l+1][0:4]
                K_temp_below = K(capital_theta(theta_below_deriv, theta_r, theta_s), K_s, m) #hydraulic conductivity based on van Genuchten parameters and theta
            except:
                K_temp_below=0

            if (l==wf_number_for_h_p):
                h_p=h_p_store
            else:
                h_p=0
            #h_p = h_p*(True if l==wf_number_for_h_p else False)
            ##h_p = h_p*(l==wf_number_for_h_p)
            #print(h_p)
    #         print('wf_number_for_h_p:')
    #         print(wf_number_for_h_p)
    #         print('l in derivs:')
    #         print(l)
    #         print(l==wf_number_for_h_p)
            #print(l)
            #print(h_p)
            #3 feb move
            ###Z, theta, layer_number, number_in_soil = current_states[l][0:4] #current states has the format of Z, theta, and then the layer number the wetting front is in
            theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[int(layer_number)] #it's useful to have the right parameters on hand; the parameters[layer_number] contains relevant parameter for layer_number. Layer number is also 0-indexed, such that the most superficial layer has a layer number of 0.
            if (l>0): #the idea here is that is l==0, then there is no current_states[l-1]. So we set Z_above, theta_above, layer_number_above to be Z, theta, layer_number in only the case of the most superficial wetting front, although calculation of derivatives for the most superficial wetting front will never need these variables.
                Z_above, theta_above, layer_number_above = self.current_states[l-1][0:3]
            else:
                Z_above, theta_above, layer_number_above = Z, theta, layer_number

            if (l<maximum_wetting_front_number): #the idea here is that the deepest wetting front, with a wetting front number of maximum_wetting_front_number, will have special conditions, as its Z value cannot change (due to the free drainage boundary condition)
                Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3] #these are the info for the wetting front below the current one, which will be useful later for calculating derivatives (namely, capillary suction calculation needs to know the moisture below the wetting front in question)
                if ((layer_number_below==layer_number)): ##layer_number_above==layer_number###the following code will set derivatives for wetting fronts which have, in the same layer, both wetting fronts above and below. Different conditions are necessary to calculate derivatives when the wetting front is either the deepest or most superficial in its layer.
                    theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[layer_number] #this is probably redundant

                    if (layer_number == 0): #currently, LGAR works with exactly 3 layers (which can be set to have the same properties). Derivatives have different forms for wetting fronts in different layers, as moisture in above layers must be considered for wetting fronts in deeper layers. Therefore dZdt calculations look similar between layers, albeit there are more terms in parts of dZdt for deeper layers. Eventually this will probably be replaced with code that generally calculates dZdt based on what layer the wetting front is in.
                        K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m) #hydraulic conductivity based on van Genuchten parameters and theta
                        psi_b = self.psi_b_vec[layer_number]
                        lambdaa = self.lambda_vec[layer_number]
                        if self.closed_form_capillary_drive_term:
                            G_temp = G_closed(theta,theta_below,psi_b,lambdaa,theta_r,theta_s)
                        else:
                            G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n,m,alpha), alpha, n, m, K_s) #capillary suction the wetting front experiences, which is in part controlled by the theta value of the wetting front below the current one, where the wetting front below will be in the same layer as the current wetting front

                        #print('G_temp closed form:')
                        #print(G_temp)
                        ###testing closed form G by commenting out line below and including 3 lines above
                        #G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n,m,alpha), alpha, n, m, K_s) #capillary suction the wetting front experiences, which is in part controlled by the theta value of the wetting front below the current one, where the wetting front below will be in the same layer as the current wetting front
                        #print('G_temp integrated:')
                        #print(G_temp)
                        #1/0
                        f_temp = [1/abs(theta-theta_below)*(K_temp+K_s*(G_temp+h_p)/Z)#, #f_temp contains the derivatives, dZdt and dthetadt, for the wetting front, and will be appended to f. again, dZdt is calculated but dthetadt is set to 0 and theta is later calculated via mass balance
                                    #(1/abs(theta-theta_below))*(K_temp+(K_temp-K_temp_below)*(G_temp+h_p)/Z)
                                  #0
                                  ]
                        f.append(f_temp)
                        #f[number_in_soil] = f_temp[0]



                    if (layer_number > 0):
                        K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m)
                        psi_b = self.psi_b_vec[layer_number]
                        lambdaa = self.lambda_vec[layer_number]
                        if self.closed_form_capillary_drive_term:
                            G_temp = G_closed(theta,theta_below,psi_b,lambdaa,theta_r,theta_s)
                        else:
                            G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n,m,alpha), alpha, n, m, K_s)

                        psi_current = psi_of_theta(theta, theta_s, theta_r, n, m, alpha) #psi is the hydraulic head. One of the central ideas of LGAR is that wetting fronts in deeper layers extend upwards to layers above, but with the same psi value (hydraulic head) rather than the same theta value. So, the same "wetting front" in terms of psi extends between multiple layers (whereas the LGAR code tracks distinct depth,theta pairs, so a wetting front going between 2 soil layers with a single psi value is represented in LGAR as 2 wetting fronts with different theta values but the same psi value). Anyway, the theta value curresponding to the wetting fornt's psi value is necessary for dZdt, so psi is calculated here and is then used to find the theta value in the layer above corresponding to this wetting front.
                        capital_theta_above_vec = []
                        K_above_vec = []
                        for k in range(1,layer_number+1):
                            capital_theta_above = capital_theta_of_psi(psi_current,alpha_vec[layer_number-k],n_vec[layer_number-k],m_vec[layer_number-k])
                            #capital_theta_above_vec.append(capital_theta_above)
                            capital_theta_above_vec = [capital_theta_above] + capital_theta_above_vec #adds capital_theta_above to start of list
                            K_temp_above = K(capital_theta_above, K_s_vec[layer_number-k],m_vec[layer_number-k])
                            #K_above_vec.append(K_temp_above)
                            K_above_vec = [K_temp_above] + K_above_vec #adds K_temp_above to start of list
                        Z_above = sum(max_depth_vec[0:layer_number])

                        #####3 jan 2023: trying new f_p and dZdt formulation
                        # sum_of_thickness_over_conductivity_terms = 0
                        # for k in range(0,(len(K_above_vec))):
                        #     sum_of_thickness_over_conductivity_terms = sum_of_thickness_over_conductivity_terms + max_depth_vec[k]/K_above_vec[k]
                        # f_temp = [(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*K_s/K_temp)*(1/(sum_of_thickness_over_conductivity_terms+Z/K_temp))#, #was 0
                        #             #(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*(K_temp-K_temp_below)/K_temp)*(1/(max_depth_vec[0]/K_temp_2_above+max_depth_vec[1]/K_temp_above+Z/K_temp))
                        #           #0
                        #           ]
                        sum_of_thickness_over_conductivity_terms = 0
                        #sum_of_thickness_over_conductivity_terms_K_s = 0
                        for k in range(0,(len(K_above_vec))):
                            sum_of_thickness_over_conductivity_terms = sum_of_thickness_over_conductivity_terms + max_depth_vec[k]/K_above_vec[k]
                            #sum_of_thickness_over_conductivity_terms_K_s = sum_of_thickness_over_conductivity_terms_K_s + max_depth_vec[k]/K_s_vec[k]
                        sum_of_thickness_over_conductivity_terms = sum_of_thickness_over_conductivity_terms + Z/K_temp
                        #sum_of_thickness_over_conductivity_terms_K_s = sum_of_thickness_over_conductivity_terms_K_s + Z/K_s
                        K_composite = (Z+Z_above)/sum_of_thickness_over_conductivity_terms
                        #K_s_composite = (Z+Z_above)/sum_of_thickness_over_conductivity_terms_K_s
                        f_temp = [K_s*G_temp/(Z+Z_above)+K_composite]

                        f.append(f_temp)


                        #####this was the old code for when LGAR had exactly 3 layers. the code above beginning with "if (layer_number > 0):" is the new code that allows for any number of layers.
                    # if (layer_number == 1): #currently, LGAR works with exactly 3 layers (which can be set to have the same properties). Derivatives have different forms for wetting fronts in different layers, as moisture in above layers must be considered for wetting fronts in deeper layers. Therefore dZdt calculations look similar between layers, albeit there are more terms in parts of dZdt for deeper layers. Eventually this will probably be replaced with code that generally calculates dZdt based on what layer the wetting front is in.
                    #     K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m)
                    #     G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n,m,alpha), alpha, n, m, K_s)
                    #     psi_current = psi_of_theta(theta, theta_s, theta_r, n, m, alpha) #psi is the hydraulic head. One of the central ideas of LGAR is that wetting fronts in deeper layers extend upwards to layers above, but with the same psi value (hydraulic head) rather than the same theta value. So, the same "wetting front" in terms of psi extends between multiple layers (whereas the LGAR code tracks distinct depth,theta pairs, so a wetting front going between 2 soil layers with a single psi value is represented in LGAR as 2 wetting fronts with different theta values but the same psi value). Anyway, the theta value curresponding to the wetting fornt's psi value is necessary for dZdt, so psi is calculated here and is then used to find the theta value in the layer above corresponding to this wetting front.
                    #     capital_theta_above = capital_theta_of_psi(psi_current,alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     K_temp_above = K(capital_theta_above, K_s_vec[layer_number-1],m_vec[layer_number-1]) #further, the hydraulic conductivity in the layer above corresponding to the current wetting front is necessary for dZdt.
                    #     Z_above = max_depth_vec[0] #the wetting front which in the layer above the current one and corresponding to the psi value for the current wetting front, by definition in LGAR, is a wetting front that spans the entire layer above and so has this depth.
                    #     f_temp = [(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*K_s/K_temp)*(1/(max_depth_vec[0]/K_temp_above+Z/K_temp))#,
                    #                 #(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*(K_temp-K_temp_below)/K_temp)*(1/(max_depth_vec[0]/K_temp_above+Z/K_temp))
                    #               #0
                    #               ]
                    #     f.append(f_temp)
                    #     #f[number_in_soil] = f_temp[0]
                    # if (layer_number == 2): #this is the same idea as with layer_number == 1, but now the info of 2 wetting fronts above the current one (corresponding to the same psi value) must be considered for dZdt
                    #     K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m)
                    #     G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n,m,alpha), alpha, n, m, K_s)
                    #     psi_current = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                    #     capital_theta_above = capital_theta_of_psi(psi_current,alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     K_temp_above = K(capital_theta_above, K_s_vec[layer_number-1],m_vec[layer_number-1])
                    #     capital_theta_2_above = capital_theta_of_psi(psi_current,alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     K_temp_2_above = K(capital_theta_above, K_s_vec[layer_number-2],m_vec[layer_number-2])
                    #     Z_above = max_depth_vec[0] + max_depth_vec[1]
                    #     f_temp = [(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*K_s/K_temp)*(1/(max_depth_vec[0]/K_temp_2_above+max_depth_vec[1]/K_temp_above+Z/K_temp))#, #was 0
                    #                 #(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*(K_temp-K_temp_below)/K_temp)*(1/(max_depth_vec[0]/K_temp_2_above+max_depth_vec[1]/K_temp_above+Z/K_temp))
                    #               #0
                    #               ]
                    #     f.append(f_temp)
                    #     #f[number_in_soil] = f_temp[0]




                elif (layer_number_below!=layer_number):#in this case, the wetting front is the deepest wetting front in its layer. This means that its dZdt value is 0 and its theta value will be updated later, with the condition that it will have the same psi value as the wetting front below it (which is in the next layer)
                    #f_temp = [0,0]
                    f_temp = [0]
                    f.append(f_temp)
                    #f[number_in_soil] = f_temp[0]

                #####
                # elif (layer_number_above!=layer_number):#In this case, the layer number of the wetting front above is not the same as the layer number of the current wetting front, meaning that the current wetting front is the most superficial in its layer. If it is not the only wetting front in its layer (checked in the next line), its derivatives are similar to the case where there is a wetting front above and below the present one, although in this case the wetting front is necessarily not in the most superficial layer. Later in the mass balance (so not in the derivs function), the psi value of the wetting front directly above this one will have its psi set to the psi value of the current wetting front.
                #     #if ( (layer_number==1) & ( sum(np.array(current_states)[:,2]==1)>=2 ) ):
                #     if ( (layer_number==1) & ( layers_vec.count(1)>=2 ) ):
                #     #if ( (layer_number==1) & ( sum(using_list_comp(current_states,[2]]))>=2 ) ):
                #         theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[int(layer_number)]
                #         K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m)
                #         G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha), alpha, n, m, K_s)
                #         current_psi_value = psi_of_theta(theta,theta_s,theta_r,n,m,alpha)
                #         theta_above = theta_of_psi(current_psi_value, theta_r_vec[0], theta_s_vec[0], alpha_vec[0], n_vec[0], m_vec[0])
                #         capital_theta_above = capital_theta(theta_above, theta_r_vec[0], theta_s_vec[0])
                #         K_temp_above = K(capital_theta_above,K_s_vec[0],m_vec[0])
                #         f_temp = [ (1/abs(theta-theta_below))*(Z+max_depth_vec[0]+(G_temp+h_p)*K_s/K_temp)*(1/(max_depth_vec[0]/K_temp_above+Z/K_temp))#,
                #                     #(1/abs(theta-theta_below))*(Z+max_depth_vec[0]+(G_temp+h_p)*(K_temp-K_temp_below)/K_temp)*(1/(max_depth_vec[0]/K_temp_above+Z/K_temp))
                #               #0
                #             ]
                #         f.append(f_temp)
                #         #f[number_in_soil] = f_temp[0]
                #     #elif ( (layer_number==1) & ( sum(np.array(current_states)[:,2]==1)==1 ) ): #in this case, while the current wetting front is the most spuerficial in its layer, it is also the only wetting front in its layer, meaning that its depth spans its whole layer and therefore dZdt = 0, and its theta value will be set such that the psi value matches that of the wetting front below this one (this psi setting will be done later, not in the derivs function)
                #     elif ( (layer_number==1) & ( layers_vec.count(1)==1 ) ):
                #         #f_temp = [0,0]
                #         f_temp = [0]
                #         f.append(f_temp)
                #         #f[number_in_soil] = f_temp[0]
                #     if (layer_number==2):#still in the case where this is the most superficial wetting front in its layer, this now just provides derivatives in the event that it is layer 2 rather than layer 1
                #         #if ( ( sum(np.array(current_states)[:,2]==2)>=2 ) ):
                #         if ( ( layers_vec.count(2)>=2 ) ):
                #             theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[layer_number]
                #             K_temp = K(capital_theta(theta, theta_r, theta_s), K_s, m)
                #             G_temp = G(psi_of_theta(theta, theta_s, theta_r, n, m, alpha), psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha), alpha, n, m, K_s)
                #             current_psi_value = psi_of_theta(theta,theta_s,theta_r,n,m,alpha)
                #             theta_above = theta_of_psi(current_psi_value, theta_r_vec[1], theta_s_vec[1], alpha_vec[1], n_vec[1], m_vec[1])
                #             capital_theta_above = capital_theta(theta_above, theta_r_vec[1], theta_s_vec[1])
                #             K_temp_above = K(capital_theta_above,K_s_vec[1],m_vec[1])
                #             theta_2_above = theta_of_psi(current_psi_value, theta_r_vec[0], theta_s_vec[0], alpha_vec[0], n_vec[0], m_vec[0])
                #             capital_theta_2_above = capital_theta(theta_2_above, theta_r_vec[0], theta_s_vec[0])
                #             K_temp_2_above = K(capital_theta_2_above,K_s_vec[0],m_vec[0])
                #             f_temp = [ (1/abs(theta-theta_below))*(Z+(max_depth_vec[0]+max_depth_vec[1])+(G_temp+h_p)*K_s/K_temp)*(1/((max_depth_vec[1])/K_temp_above+max_depth_vec[0]/K_temp_2_above+Z/K_temp))#, #0, #this should be based on new_equation
                #                         #(1/abs(theta-theta_below))*(Z+(max_depth_vec[0]+max_depth_vec[1])+(G_temp+h_p)*(K_temp-K_temp_below)/K_temp)*(1/((max_depth_vec[1])/K_temp_above+max_depth_vec[0]/K_temp_2_above+Z/K_temp))
                #                   #0
                #                 ]
                #             f.append(f_temp)
                #             #f[number_in_soil] = f_temp[0]
                #####


            elif (l==maximum_wetting_front_number):#finally, these are the derivs for the deepest wetting front. While here a nonzero deivative for Z is given, later what will happen is this derivative will be used to calculate the free drainage demand, which is subtracted from the most superficial wetting front, and the Z value for the deepest wetting front remains unchanged.
                theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[int(layer_number)]
                ###free_drainage_factor was an old idea we implemented before we were using ET; it is a way to get bottom fluxes, when nonzero, to match well between LGAR and HYDRUS. My strong recommendation is to not use this approach and instead rely on large ET values to make LGAR and HYDRUS compare well.
                #f_temp = [free_drainage_factor*K(capital_theta(theta, theta_r, theta_s),K_s,m)/theta,0]
                f_temp = [0]
                if (h_p>0): #added the 0* in 0*K_s*h_p/(sum(max_depth_vec)) on 15 June 2022
                    #f_temp = [ 0*(K(capital_theta(theta, theta_r, theta_s),K_s,m) + 0*K_s*h_p/(sum(max_depth_vec)) ) /theta,0] ### 10 dec the *0 is to get h_p to not affect free drainage, currently not there but can make h_p*0 to do this ###
                    f_temp = [ 0 ]
                ###the following line makes the lower BC not free drainage but constant flux instead
                #f_temp = [free_drainage_demand_constant/time_step,0]
                f.append(f_temp)
                #f[number_in_soil] = f_temp[0]
        #print(f)
        return(f)

    def theta_of_psi(self, psi, theta_r, theta_s, alpha, n, m): #gives volumetric water content (unitless) based on hydraulic head
        return(  theta_r + (theta_s-theta_r)/(1+(alpha*psi)**n)**m  )

     #used to calculate wall clock runtime for computations


    #for i in (range(0,1)):
    def run_model(self,length_of_simulation):

        derivs = self.derivs

        precip_data = self.precip_data
        PET_data = self.PET_data

        config = self.config
        params = config.params
        parameters = params.parameters
        self.closed_form_capillary_drive_term = config.closed_form_capillary_drive_term

        theta_r_vec   = params.theta_r_vec
        theta_s_vec   = params.theta_s_vec
        K_s_vec       = params.K_s_vec
        alpha_vec     = params.alpha_vec
        n_vec         = params.n_vec
        m_vec         = params.m_vec
        max_depth_vec = params.max_depth_vec

        relative_moisture_at_which_PET_equals_AET = params.relative_moisture_at_which_PET_equals_AET
        boundary_depths = params.boundary_depths
        max_layer = params.max_layer
        time_steps_to_record_profile = config.time_steps_to_record_profile

        time_step = config.time_step

        self.lambda_vec = []
        self.psi_b_vec = []

        for u in range(0,len(theta_r_vec)):
            p=1+2/m_vec[u]
            self.lambda_vec.append(2/(p-3))
            self.psi_b_vec.append( (p+3)*(147.8+8.1*p+0.092*p**2) / (2*alpha_vec[u]*p*(p-1)*(55.6+7.4*p+p**2)) )



        self.runtime_0 = datetime.now()
        for i in (range(self.current_time_step,length_of_simulation)):
            try:
                self.precip_data[i]
            except:
                self.precip_data = np.array(self.forcing_data['P(mm/h)'])[0:(self.current_time_step+1)] #these could probably happen in init too
                self.PET_data = np.array(self.forcing_data['PET(mm/h)'])[0:(self.current_time_step+1)]
                precip_data = self.precip_data
                PET_data = self.PET_data

            #self.previous_states = np.array(self.current_states)
            #previous_states = list(current_states)[0:len(current_states)]
            #previous_states = copy.deepcopy(current_states)
            #previous_states = deque(current_states)
            #previous_states = list(current_states)
            #previous_states = current_states

            ###seems to work-ish and is the fastest option. One of the biggest issues for speed is how previous states is set; above a lot of different options were tried and the one below seems to be fastest.
            empty_list = []
            for jj in range(0,len(self.current_states)):
                empty_list.append([self.current_states[jj][0],self.current_states[jj][1],self.current_states[jj][2],jj,i])
            self.previous_states = empty_list
            ###


            ###seems to work-ish
            # empty_array = np.empty((len(current_states),5))
            # for jj in range(0,len(current_states)):
            #     empty_array[jj]=[current_states[jj][0],current_states[jj][1],current_states[jj][2],jj,i]
            # #print(empty_list)
            # previous_states = empty_array
            ###

##### moving to the base BMI class, because this is finalization
            # if (self.first_time_step==False):
            #     if ( (i%mod_number)==0 ):
            #         ###this code saves results to the disk; it is more efficient to do this in batches rather than all at once, especially when the run is long.
            #         #print('len')
            #         #print(len(forcing_data[(i-100):i]))
            #
            #
            #
            #         # for i in (range(0,len(states_array))):
            #         #     mass_bal_vec = np.append(mass_bal_vec, (calc_mass_bal((states_array[states_array[:,4]==i]))))
            #         # mass_balance_error_vec = mass_bal_vec[0:length_of_simulation]-mass_bal_vec[0]+fluxes_cumulative[0:length_of_simulation]+h_p_vec[0:length_of_simulation]
            #
            #         forcing_data_save = forcing_data[(i-mod_number):i]
            #         forcing_data_save['runoff[mm/h]'] = self.runoff_vec
            #         forcing_data_save['actual_infil[mm/h]'] = self.actual_infiltration_vec
            #         #h_p_vec = h_p_vec[:len(precip_data)]
            #         forcing_data_save['ponded_head[mm]'] = self.h_p_vec
            #         forcing_data_save['bottom_flux[mm/h]'] = self.free_drainage_flux_vec
            #         #forcing_data_save['bottom_demand[mm]'] = free_drainage_vec
            #         forcing_data_save['water_in_soil[mm]'] = self.mass_vec
            #         forcing_data_save['mass_bal_error(mm)'] = self.mass_bal_vec#mass_balance_error_vec
            #         forcing_data_save['actual_ET_per_step(mm)'] = self.actual_ET_vec
            #
            #
            #         file_str = "pickle_jar/LGAR_output"+str(self.pickle_step)+".pkl"
            #         forcing_data_save.to_pickle(file_str)
            #
            #         # np.save('states_array',states_array)
            #         #
            #         # states_array = np.array(current_states)
            #         # total_num_wetting_fronts = len(current_states)
            #         # wetting_front_numbers = np.array(range(0,total_num_wetting_fronts))
            #         #
            #         # states_array = np.column_stack((np.array(current_states), wetting_front_numbers))#adds wetting front number; code is set up in such a way that each wetting front gets an integer number, starting at 0, from superficial to deep, just to represent the total number of wetting fronts
            #         #
            #         # states_array = np.column_stack((states_array, np.zeros(total_num_wetting_fronts)))#adds time step of 0; when current_states is updated later, the time step will be added instead of 0
            #
            #         self.pickle_step=self.pickle_step+1
            #
            #
            #
            #         #pickle.dump(forcing_data_save,open( "LGAR_output.pkl", "ab" ))
            #
            #
            #
            #         self.fluxes=[]
            #         #fluxes_cumulative=[]
            #         #f_p_vec=[]
            #         self.runoff_vec=[]
            #         self.actual_infiltration_vec=[]
            #         self.h_p_vec=[]
            #         self.free_drainage_flux_vec=[]
            #         #free_drainage_vec=[]
            #         self.mass_vec=[]
            #         self.actual_ET_vec=[]
            #         self.mass_bal_vec=[]
#####


                # print('np array of current_states:')
                # print(np.array(self.current_states))
                # print('looking at depths of current_states:')
                # print(np.array(self.current_states)[:,0])
                # print('looking at thetas of current_states:')
                # print(np.array(self.current_states)[:,1])
                # print(' ')

                #print(np.array(self.current_states).itemsize)
                # print(type(precip_mass_to_add))
                # print(precip_mass_to_add)

            #print(wilting_point_vec)

            #print(boundary_depths)
            ok_to_increment_h_p = 1
            #determines wetting front from which to subtract free drainage demand. This is also the wf that experiences infiltration from precip and also ponding head.
            wf_that_supplies_free_drainage_demand = 0
            layers_vec = using_list_comp(self.current_states,[2])[0]
            for h in range(0,len(self.parameters)):
                if (layers_vec.count(h)!=1):
                    break
                else:
                    wf_that_supplies_free_drainage_demand = wf_that_supplies_free_drainage_demand + 1
            if (wf_that_supplies_free_drainage_demand>max(self.wetting_front_numbers)):
                wf_that_supplies_free_drainage_demand = wf_that_supplies_free_drainage_demand - 1
            ###conveniently, the wf_that_supplies_free_drainage_demand is also the wf that should get fed by precip; the most superfical (in the context of LGAR this also means the wetting front with the psi value that is closest to 0) wetting front is the one that satisfies the free drainage demand because water is the least tightly held to the soil in this wetting front


            #2 feb takeout
            self.mass_vec = np.append(self.mass_vec, 0*self.h_p+(calc_mass_bal((self.previous_states)))) #this calculates the total amount of water in the soil profile for the most recent entry in states array, which will be used later to correct a potential mass balance error that would occur if calculated infiltration is enough to cause theta>theta_s for a given layer
            #15 march added 0*

            #print(h_p)
            #print('wf_that_supplies_free_drainage_demand:')
            #print(wf_that_supplies_free_drainage_demand)
            temp_derivs = derivs(self.current_states, self.h_p, wf_that_supplies_free_drainage_demand) #temp_derivs is the array containing the time derivatives for Z for each wetting front
            self.number_of_wetting_fronts = len(self.current_states) #this is the number of wetting fronts currently in the soil profile
            maximum_wetting_front_number = self.number_of_wetting_fronts-1 #because the wetting fronts are 0-indexed, this is the number of the deepest wetting front
            #wetting_front_numbers = np.array(range(0,number_of_wetting_fronts)) #this is an array containing all wetting front numbers, starting at 0, where wetting front 0 is the most superficial one

            self.wetting_front_numbers = np.arange(0,self.number_of_wetting_fronts)
            layer_number_below=len(parameters)-1 #later, the fluxes for wetting fronts will be calculated from the bottom up (after the line: for l in range(number_of_wetting_fronts-1,-1,-1): ), and so the layer below the current wetting front will begin as the deepest possible layer, with special derivatives and mass balance calculations for the deepest wetting front

            #current_mass_in_soil = calc_mass_bal(current_states)



            precip_mass_to_add = 0 #precip_mass_to_add is the amount of precipitation that enters the soil; if it is necessary later a nonzero value will be calculated
            runoff = 0 #the amount of water that leaves the system as runoff

            #here the free drainage is calculated. The free drainage demand will ultimately be subtracted from the wetting front which has the psi value closest to 0. See the last elif statement in the derivs function for more details about how this works
            theta_value_for_free_drainage = self.current_states[-1][1]
            free_drainage_demand = temp_derivs[-1][0]*config.time_step*theta_value_for_free_drainage
            #free_drainage_demand = temp_derivs[-1]*time_step*theta_value_for_free_drainage
            ###un comment for constant flux
            #free_drainage_demand = free_drainage_demand_constant

            #free_drainage_augmented = 0
            #free_drainage_vec.append(free_drainage_demand)

            #here, f_p, or the potential infiltration capacity, will be calculated in the event that both precip in the current and previous time step were both greater than 0.
            #this is due to the fact that new wetting front creation, occurring when the precip for the current time step is nonzero and the precip for the previous time step was 0, is handled separately. It might be wise to eventually change it so that a new wettinf front gets created any time there is a change in precip (rather than any time precip starts), but for now this seems to do well ...
            #this does mean that at present, the code will not simulate a new wetting front if the precip at time 0 is nonzero, and instead precip at time zero adds to the currently existing most superficial wetting front.

            ###aslo note that this code presently only works with h_p=0 and h_p_max = 0. Or at least, h_p_max was developed back before ET was added, and surface ponding has not yet been tested with ET.

            #note that the calculation of f_p is somewhat, but not very, computationally expensive; this is because it uses the function G. Any time this function is called, it performs trapezoid rule integration, which is a bit expensive. Therefore, f_p is only calculated when it is relevant, that is if there is nonzero precip in both the current and previous time steps.
            if (self.h_p==0):
                if (i>0):
                    if ( (self.precip_data[i]>0)&(self.precip_data[i-1]>0) ): #necessary condition for the calculation of f_p, implying that there is not a new wetting front
                        wf_that_supplies_free_drainage_demand_data = self.current_states[wf_that_supplies_free_drainage_demand]
                        Z_fp, theta_fp, layer_fp = wf_that_supplies_free_drainage_demand_data[0:3]
                        K_s_temp = self.K_s_vec[layer_fp]
                        if (len(self.current_states)==len(self.parameters)): #this is in the case of free drainage, where there is no capillary suction
                            G_fp=0
                        else:
                            wf_below_free_drainage_one = self.current_states[wf_that_supplies_free_drainage_demand+1]
                            Z_below_fp, theta_below_fp, layer_below_fp = wf_below_free_drainage_one[0:3]
                            psi_b = self.psi_b_vec[layer_fp]
                            lambdaa = self.lambda_vec[layer_fp]
                            theta_r = self.theta_r_vec[layer_fp]
                            theta_s = self.theta_s_vec[layer_fp]
                            if self.closed_form_capillary_drive_term:
                                G_fp = G_closed(theta_s,theta_below_fp,psi_b,lambdaa,theta_r,theta_s)
                            else:
                                G_fp = G(psi_of_theta(self.theta_s_vec[layer_fp], self.theta_s_vec[layer_fp], self.theta_r_vec[layer_fp], self.n_vec[layer_fp], self.m_vec[layer_fp], self.alpha_vec[layer_fp]), psi_of_theta(theta_below_fp, self.theta_s_vec[layer_fp], self.theta_r_vec[layer_fp], self.n_vec[layer_fp], self.m_vec[layer_fp], alpha_vec[layer_fp]), self.alpha_vec[layer_fp], self.n_vec[layer_fp], self.m_vec[layer_fp], self.K_s_vec[layer_fp])


                        # K_temp_below=0
                        # try:
                        #     Z_below_deriv, theta_below_deriv, layer_number_below_deriv, number_in_soil_below_deriv = self.current_states[wf_that_supplies_free_drainage_demand+1][0:4]
                        #     if (layer_number_below_deriv==layer_fp):
                        #         K_temp_below = K(capital_theta(theta_below_deriv, theta_r_vec[layer_fp], theta_s_vec[layer_fp]), K_s_vec[layer_fp], m_vec[layer_fp]) #hydraulic conductivity based on van Genuchten parameters and theta
                        # except:
                        #     K_temp_below=0
                        #
                        # K_temp = K(capital_theta(theta_fp, theta_r_vec[layer_fp], theta_s_vec[layer_fp]), K_s_vec[layer_fp], m_vec[layer_fp])


                        if (layer_fp==0): #the calculation of potential infiltration capacity has different forms for different layers, closely analgous to the different forms for dZdt for the different layers in the derivs function.
                            f_p = K_s_temp*(1+G_fp/Z_fp)



                        if (layer_fp > 0):
                            sum_of_thickness_over_saturated_conductivity_terms = 0
                            Z_above_for_fp = 0
                            for k in range(0,layer_fp):
                                Z_above_for_fp = Z_above_for_fp + self.max_depth_vec[k]
                                sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + self.max_depth_vec[k]/self.K_s_vec[k]

                            #####3 jan 2023: trying new f_p and dZdt formulation
                            # f_p = (Z_fp+Z_above_for_fp+G_fp)/(sum_of_thickness_over_saturated_conductivity_terms+Z_fp/self.K_s_vec[layer_fp])
                            # if (layer_fp==1):
                            #     f_p = (Z_fp+self.max_depth_vec[0]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+Z_fp/self.K_s_vec[1])
                            # if (layer_fp==2):
                            #     f_p = (Z_fp+self.max_depth_vec[0]+self.max_depth_vec[1]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+self.max_depth_vec[1]/self.K_s_vec[1]+Z_fp/self.K_s_vec[2])

                            sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + Z_fp/K_s_temp
                            K_s_composite = (Z_fp+Z_above_for_fp)/(sum_of_thickness_over_saturated_conductivity_terms)
                            f_p = K_s_composite+K_s_temp*G_fp/(Z_fp+Z_above_for_fp)
                            #f_p = K_s_composite+K_s_composite*G_fp/(Z_fp+Z_above_for_fp)
                            #f_p = K_s_temp*(1+G_fp/(Z_fp+Z_above_for_fp))

                            if ((layer_fp==(len(self.parameters)-1))&(self.current_states[0][1]==self.theta_s_vec[0])&(len(self.current_states)==len(self.parameters))): #
                                f_p = 0#self.K_s_vec[-1]# that was for free drainage, edited 23 march

                        # if (layer_fp==1):
                        #     f_p = (Z_fp+self.max_depth_vec[0]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+Z_fp/self.K_s_vec[1])
                        # if (layer_fp==2):
                        #     f_p = (Z_fp+self.max_depth_vec[0]+self.max_depth_vec[1]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+self.max_depth_vec[1]/self.K_s_vec[1]+Z_fp/self.K_s_vec[2])
                        #
                        # if ((layer_fp==2)&(self.current_states[0][1]==self.theta_s_vec[0])&(len(self.current_states)==len(parameters))): #
                        #     f_p = 0#self.K_s_vec[-1]# that was for free drainage, edited 23 march

                        #f_p_vec.append(f_p)

                        if (precip_data[i]<(f_p+free_drainage_demand/time_step) ): #this code ensures that the correct amount of water infiltrates to the soil: f_p*time_step + free drainage demand if precip data is sufficiently larger, and just precip_data[i]*time_step otherwise
                            precip_mass_to_add = (precip_data[i])*time_step
                            if ( (self.current_states[0][1]==theta_s_vec[0]) & (len(parameters)==len(self.current_states)) ): ###trying on 29 dec ###
                                precip_mass_to_add = free_drainage_demand
                                runoff = precip_data[i]*time_step - precip_mass_to_add
                        else:
                            precip_mass_to_add = (f_p)*time_step+free_drainage_demand
                            ###14 september 2022: testing creation of new wetting front if precip>f_p
                            # if (self.current_states[0][1] < theta_s_vec[0]):
                            #     precip_mass_to_add_new_wf_due_to_intense_precip = precip_mass_to_add
                            #     #precip_mass_to_add = 0
                            ###
                            if ((layer_fp==(len(parameters)-1))&(len(self.current_states)==len(parameters))): #June 16 -- if mass balance errors occur, look here first
                                precip_mass_to_add = precip_mass_to_add #+ free_drainage_demand ###9 September addition to allow for complete vadose zone saturation in intense precip; added + precip_mass_to_add
                            runoff = precip_data[i]*time_step - precip_mass_to_add ###30 nov### runoff + ### 10 dec extra indent

                else: #in this special case, i==0. In the event that there is precip at time 0, infiltration is added to the currently existing most superficial wetting front, rather than a new wetting front being created. Therefore, f_p is calculated, whether it is needed or not, for just the first time step.
                    wf_that_supplies_free_drainage_demand_data = self.current_states[wf_that_supplies_free_drainage_demand]
                    Z_fp, theta_fp, layer_fp = wf_that_supplies_free_drainage_demand_data[0:3]
                    K_s_temp = params.K_s_vec[int(layer_fp)]
                    if (len(self.current_states)==len(parameters)): #this is in the case of free drainage, where there is no capillary suction
                        G_fp=0
                    else:
                        wf_below_free_drainage_one = self.current_states[wf_that_supplies_free_drainage_demand+1]
                        Z_below_fp, theta_below_fp, layer_below_fp = wf_below_free_drainage_one[0:3]
                        psi_b = self.psi_b_vec[layer_fp]
                        lambdaa = self.lambda_vec[layer_fp]
                        theta_r = self.theta_r_vec[layer_fp]
                        theta_s = self.theta_s_vec[layer_fp]
                        if self.closed_form_capillary_drive_term:
                            G_fp = G_closed(theta_s,theta_below_fp,psi_b,lambdaa,theta_r,theta_s)
                        else:
                            G_fp = G(psi_of_theta(theta_s_vec[layer_fp], theta_s_vec[layer_fp], theta_r_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], alpha_vec[layer_fp]), psi_of_theta(theta_below_fp, theta_s_vec[layer_fp], theta_r_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], alpha_vec[layer_fp]), alpha_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], K_s_vec[layer_fp])

                    if (layer_fp==0): #the calculation of potential infiltration capacity has different forms for different layers, closely analgous to the different forms for dZdt for the different layers in the derivs function.
                        f_p = K_s_temp*(1+G_fp/Z_fp)

                    if (layer_fp > 0):
                        sum_of_thickness_over_saturated_conductivity_terms = 0
                        Z_above_for_fp = 0
                        for k in range(0,layer_fp):
                            Z_above_for_fp = Z_above_for_fp + self.max_depth_vec[k]
                            sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + self.max_depth_vec[k]/self.K_s_vec[k]

                        #####3 jan 2023: trying new f_p and dZdt formulation
                        #f_p = (Z_fp+Z_above_for_fp+G_fp)/(sum_of_thickness_over_saturated_conductivity_terms+Z_fp/self.K_s_vec[layer_fp])

                        sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + Z_fp/K_s_temp
                        K_s_composite = (Z_fp+Z_above_for_fp)/(sum_of_thickness_over_saturated_conductivity_terms)
                        f_p = K_s_composite+K_s_temp*G_fp/(Z_fp+Z_above_for_fp)
                        #f_p = K_s_composite+K_s_composite*G_fp/(Z_fp+Z_above_for_fp)
                        #f_p = K_s_temp*(1+G_fp/(Z_fp+Z_above_for_fp))

                        # if (layer_fp==1):
                        #     f_p = (Z_fp+self.max_depth_vec[0]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+Z_fp/self.K_s_vec[1])
                        # if (layer_fp==2):
                        #     f_p = (Z_fp+self.max_depth_vec[0]+self.max_depth_vec[1]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+self.max_depth_vec[1]/self.K_s_vec[1]+Z_fp/self.K_s_vec[2])

                        if ((layer_fp==len(parameters))&(self.current_states[0][1]==self.theta_s_vec[0])&(len(self.current_states)==len(parameters))): #
                            f_p = 0#self.K_s_vec[-1]# that was for free drainage, edited 23 march

                    ###experimental 15 march takeout. this was for free drainage.
                    # if ((layer_fp==2)&(len(self.current_states)==len(parameters))): #in this case, there is only one psi value for all wetting fronts and so free drainage dictates f_p. It might be a wise idea to later change it so that a new wetting front gets generated whenever there is a change in precip intensity rather than just a change in precip ...
                    #     f_p = params.K_s_vec[-1]#*current_states[-1][1] ### commented out theta component 29 nov


                    # if ((layer_fp==2)&(self.current_states[0][1]==self.theta_s_vec[0])&(len(self.current_states)==len(self.parameters))):
                    #     f_p = 0



                    #f_p_vec.append(f_p)

                    if (precip_data[i]<(f_p+free_drainage_demand/time_step) ): #this code ensures that the correct amount of water infiltrates to the soil: f_p + free drainage demand if precip data is sufficiently larger, and just precip_data[i]*time_step otherwise
                        precip_mass_to_add = (precip_data[i])*time_step
                    else:
                        precip_mass_to_add = (f_p)*time_step+free_drainage_demand
                        runoff = precip_data[i]*time_step - precip_mass_to_add ###30 nov### runoff +

                    if (precip_data[i]>0)&(precip_data[i-1]==0):#eh good to include just to be safe but shouldn't be necessary b/c of conditions for making new wetting front
                        precip_mass_to_add = 0

                #this is probably redundant as well, just makes sure that the precip_mass_to_add is not physically impossible
                #but removed on 23 march because I think it might be causing an error, because this check does not consider that precip_mass_to_add could be greater than precip*time_step if h_p>0
                # if (precip_mass_to_add>(precip_data[i]*time_step)):
                #     precip_mass_to_add = precip_data[i]*time_step
                #     runoff = 0


            if ( self.h_p>0 ):
                wf_that_supplies_free_drainage_demand_data = self.current_states[wf_that_supplies_free_drainage_demand]
                Z_fp, theta_fp, layer_fp = wf_that_supplies_free_drainage_demand_data[0:3]
                K_s_temp = self.K_s_vec[layer_fp]
                if (len(self.current_states)==len(self.parameters)): #this is in the case of free drainage, where there is no capillary suction
                    G_fp=0
                else:
                    wf_below_free_drainage_one = self.current_states[wf_that_supplies_free_drainage_demand+1]
                    Z_below_fp, theta_below_fp, layer_below_fp = wf_below_free_drainage_one[0:3]
                    psi_b = self.psi_b_vec[layer_fp]
                    lambdaa = self.lambda_vec[layer_fp]
                    theta_r = self.theta_r_vec[layer_fp]
                    theta_s = self.theta_s_vec[layer_fp]
                    if self.closed_form_capillary_drive_term:
                        G_fp = G_closed(theta_s,theta_below_fp,psi_b,lambdaa,theta_r,theta_s)
                    else:
                        G_fp = G(psi_of_theta(theta_s_vec[layer_fp], theta_s_vec[layer_fp], theta_r_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], alpha_vec[layer_fp]), psi_of_theta(theta_below_fp, theta_s_vec[layer_fp], theta_r_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], alpha_vec[layer_fp]), alpha_vec[layer_fp], n_vec[layer_fp], m_vec[layer_fp], K_s_vec[layer_fp])

                if (layer_fp==0): #the calculation of potential infiltration capacity has different forms for different layers, closely analgous to the different forms for dZdt for the different layers in the derivs function.
                    f_p = K_s_temp*(1+(G_fp+self.h_p)/Z_fp)

                if (layer_fp > 0):
                    sum_of_thickness_over_saturated_conductivity_terms = 0
                    Z_above_for_fp = 0
                    for k in range(0,layer_fp):
                        Z_above_for_fp = Z_above_for_fp + self.max_depth_vec[k]
                        sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + self.max_depth_vec[k]/self.K_s_vec[k]

                    #####3 jan 2023: trying new f_p and dZdt formulation
                    #f_p = (Z_fp+Z_above_for_fp+G_fp+self.h_p)/(sum_of_thickness_over_saturated_conductivity_terms+Z_fp/self.K_s_vec[layer_fp])

                    sum_of_thickness_over_saturated_conductivity_terms = sum_of_thickness_over_saturated_conductivity_terms + Z_fp/K_s_temp
                    K_s_composite = (Z_fp+Z_above_for_fp)/(sum_of_thickness_over_saturated_conductivity_terms)
                    #f_p = K_s_composite+K_s_temp*G_fp/(Z_fp+Z_above_for_fp)
                    #f_p = K_s_composite+K_s_composite*G_fp/(Z_fp+Z_above_for_fp)
                    f_p = K_s_temp*(1+G_fp/(Z_fp+Z_above_for_fp))

                    # if (layer_fp==1):
                    #     f_p = (Z_fp+self.max_depth_vec[0]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+Z_fp/self.K_s_vec[1])
                    # if (layer_fp==2):
                    #     f_p = (Z_fp+self.max_depth_vec[0]+self.max_depth_vec[1]+G_fp)/(self.max_depth_vec[0]/self.K_s_vec[0]+self.max_depth_vec[1]/self.K_s_vec[1]+Z_fp/self.K_s_vec[2])

                if ((layer_fp==(len(parameters)-1))&(self.current_states[0][1]==self.theta_s_vec[0])&(len(self.current_states)==len(parameters))): #
                    f_p = 0#self.K_s_vec[-1]# that was for free drainage, edited 23 march

                hp_temp = self.h_p
                if ((layer_fp==(len(self.parameters) -1) )&(len(self.current_states)==len(self.parameters))):
                    self.h_p = self.h_p + self.precip_data[i]*time_step - f_p*time_step - free_drainage_demand*0 #15 march, - f_p*time_step was - f_p*time_step*0
                else:
                    self.h_p = self.h_p + self.precip_data[i]*time_step - f_p*time_step - free_drainage_demand*0

                #print(self.h_p)
                if (self.h_p>self.h_p_max):
                    runoff = runoff + (self.h_p - self.h_p_max)
                    self.h_p = self.h_p_max

                if (self.h_p<0):
                    precip_mass_to_add = hp_temp + self.precip_data[i]*time_step
                    self.h_p = 0
                    ok_to_increment_h_p = 0

                if (self.h_p>0):
                    if ((layer_fp==(len(self.parameters)-1) )&(len(self.current_states)==len(self.parameters))):
                        precip_mass_to_add = f_p*time_step + free_drainage_demand*0
                    else:
                        precip_mass_to_add = f_p*time_step + free_drainage_demand*0

                    # if (self.current_states[0][1]==self.theta_s_vec[0]):
                    #     precip_mass_to_add = 0

                # print(' ')
                # print('f_p used for h_p calc')
                # print(f_p)
                # print('new h_p')
                # print(self.h_p)
                # print(' ')


            ###here the ET is calculated
            #old Bydyko method
            #actual_ET_demand = PET_data[i]*((self.current_states[0][1]-self.wilting_point_vec[0])/(relative_moisture_at_which_PET_equals_AET*theta_s_vec[0]-self.wilting_point_vec[0]))*time_step

            #new van Genuchten method
            #1/(1+(h/h_50)^3)
            theta_fc = (self.theta_s_vec[0]-self.theta_r_vec[0])*relative_moisture_at_which_PET_equals_AET+self.theta_r_vec[0]
            #(theta,theta_s,theta_r,n,m,alpha)
            #h_fc = psi_of_theta(theta_fc,self.theta_s_vec[0],self.theta_r_vec[0],self.n_vec[0],self.m_vec[0],self.alpha_vec[0])
            theta_50 = (theta_fc-self.wilting_point_vec[0])*1/2+self.wilting_point_vec[0]
            h_50 = psi_of_theta(theta_50,self.theta_s_vec[0],self.theta_r_vec[0],self.n_vec[0],self.m_vec[0],self.alpha_vec[0])
            h = psi_of_theta(self.current_states[0][1],self.theta_s_vec[0],self.theta_r_vec[0],self.n_vec[0],self.m_vec[0],self.alpha_vec[0])

            actual_ET_demand = PET_data[i]*(1/(1+(h/h_50)**3))*time_step



            # if (self.current_states[0][2]!=self.current_states[1][2]):
            #     theta_fc_1 = (self.theta_s_vec[1]-self.theta_r_vec[1])*relative_moisture_at_which_PET_equals_AET+self.theta_r_vec[1]
            #     #(theta,theta_s,theta_r,n,m,alpha)
            #     #h_fc = psi_of_theta(theta_fc,self.theta_s_vec[0],self.theta_r_vec[0],self.n_vec[0],self.m_vec[0],self.alpha_vec[0])
            #     theta_50_1 = (theta_fc-self.wilting_point_vec[1])*1/2+self.wilting_point_vec[1]
            #     h_50_1 = psi_of_theta(theta_50_1,self.theta_s_vec[1],self.theta_r_vec[1],self.n_vec[1],self.m_vec[1],self.alpha_vec[1])
            #     h_50_1 = h_50
            #     h_1 = psi_of_theta(self.current_states[1][1],self.theta_s_vec[1],self.theta_r_vec[1],self.n_vec[1],self.m_vec[1],self.alpha_vec[1])
            #     frac_in_1 = self.current_states[1][0]/(self.current_states[1][1]+self.max_depth_vec[0])
            #     frac_in_0 = self.max_depth_vec[0]/(self.current_states[1][1]+self.max_depth_vec[0])
            #     actual_ET_demand = PET_data[i]*(1/(1+(h/h_50)**3))*(frac_in_0)*time_step + PET_data[i]*(1/(1+(h_1/h_50_1)**3))*(frac_in_1)*time_step






            if (actual_ET_demand<0):
                actual_ET_demand=0
            if ((actual_ET_demand)>(PET_data[i]*time_step)):
                actual_ET_demand = PET_data[i]*time_step

            #print('actual ET demand (mm)')
            #print(actual_ET_demand)

            ###17 March: adding ET from h_p, in such a way that h_p supplies the part of PET that the soil can't
            # h_p_temp = self.h_p
            # actual_ET_demand_from_h_p = 0
            # if ( (actual_ET_demand<(PET_data[i]*time_step)) & (self.h_p>0) ):
            #     actual_ET_demand_from_h_p = PET_data[i]*time_step - actual_ET_demand
            #     self.h_p = self.h_p - actual_ET_demand_from_h_p
            #     if (self.h_p<0):
            #         self.h_p = 0
            #         actual_ET_demand_from_h_p = h_p_temp

            # if (actual_ET_demand_from_h_p>0):
            #     print('actual_ET_demand_from_h_p')
            #     print(actual_ET_demand_from_h_p)
            # if (actual_ET_demand_from_h_p>0):
            #     1/0

            ###this is a loop that updates Z and theta for each wetting front. It updates these starting with the deepest wetting front and moving towards the most superficial.
            for l in range(self.number_of_wetting_fronts-1,-1,-1):

                time_step = config.time_step

                theta_r_vec=self.theta_r_vec #yeah for example don't need these, not super clean
                theta_s_vec   = params.theta_s_vec
                K_s_vec       = params.K_s_vec
                alpha_vec     = params.alpha_vec
                n_vec         = params.n_vec
                m_vec         = params.m_vec
                max_depth_vec = params.max_depth_vec


                #current_states[l][0] = current_states[l][0] + temp_derivs[l][0]*time_step #I left this in as a comment because it shows in generaly how the folling code works; current_states[l][0] is the Z value for wetting front l and is made greater via the derivative for Z for that wetting front. This line will be implemented in the following conditionals when necessary (it will not be necessary for the deepest wetting front)
                layer_number = self.current_states[l][2]
                layer_number = int(layer_number)
                theta_r, theta_s, K_s, alpha, n, m, max_depth = parameters[int(layer_number)]

                if (l==max(self.wetting_front_numbers)): #this means the deepest wetting front
                    theta_value_for_free_drainage = self.current_states[l][1]
                    free_drainage_demand = temp_derivs[l][0]*time_step*theta_value_for_free_drainage #see the last elif statement in the derivs function for more info on how free drainage works
                    #free_drainage_demand = temp_derivs[l]*time_step*theta_value_for_free_drainage
                    ###un comment for constant flux
                    #free_drainage_demand = free_drainage_demand_constant

                Z, theta, layer_number = self.current_states[l][0:3]
                layer_number = int(layer_number)
                if (l>0):
                    layer_number_above = self.current_states[l-1][2]
                else: #for the case of the most superficial wetting front, there is no soil layer above, and so we artifically set layer_number_above = layer_number if l==0, meaning the most superficial wetting front. This is just to prevent potential errors; in reality, the most superficial wetting front needs no "layer above" for its derivative calculation.
                    layer_number_above = layer_number
                if (l<maximum_wetting_front_number):
                    layer_number_below = self.current_states[l+1][2] #similarly, this code determines the layer number of the wetting front that is below the current wetting front, with a special case for the deepest wetting front
                if (l==maximum_wetting_front_number):
                    layer_number_below = layer_number + 1

                ###layer_number is the layer number that the current wetting front is in, where the current layer is represented by current_states[l]. layer_number_above is the layer number that the wetting front above the current one is in. This generally could either be the same layer as current_layer or the one above. layer_number_below is the layer number of the wetting front that is below the current one. generally if can be the same as the current wetting front layer or one greater. a wetting front is the deepest wetting front in its layer if layer_number_below!=layer_number.
                ###in the following section, the theta value for each wetting front is updated based on mass balance and dZdt values.
                if ( (l<maximum_wetting_front_number)&(layer_number_below==layer_number) ): #now begins the long section of the code where theta values are calculated via mass balance. In this section, wetting fronts that have wetting fronts both above and below them in the same layer, and the most superficial wetting front, will have their theta values calculated.
                    if ((layer_number==0)): #in the case of wetting fronts in the top layer, the mass balance is fairly straightforeard.
                        #previous_states = states_array[states_array[:,4]==(len(states_array)-1)] #first, the previous_states, or the states from the end of the previous time step, are made as a variable.
                        #previous_states = states_array[states_array[:,4]==i]
                        #previous_states = states_array[states_array[:,4]==(len(states_array)-1)]
                        prior_mass = self.previous_states[l][0]*(self.previous_states[l][1]-self.previous_states[l+1][1]) - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0) #then, prior_mass is calculated. This is actually the quantity of the water in the wetting front from the previous time step, but also with relevant fluxes (free draiange and infiltration) added in. Therefore, it must be the case that the mass of water in this wetting front must be equal to the mass of water in the wetting front previously plus any relevant changes.

                        self.current_states[l][0] = self.current_states[l][0] + temp_derivs[l][0]*time_step #new Z value is recorded
                        #current_states[l][0] = current_states[l][0] + temp_derivs[l]*time_step
                        self.current_states[l][1] = min(theta_s_vec[layer_number], prior_mass/self.current_states[l][0]+self.current_states[l+1][1]) #finally, new theta value is recorded such that the wetting front conserves mass. Again, this is considering relevant fluxes.
                        Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]#might be necessary for BMI update implementation







                    if (layer_number > 0):
                        self.current_states[l][0] = self.current_states[l][0] + temp_derivs[l][0]*time_step

                        Z, theta, layer_number = self.current_states[l][0:3]
                        layer_number = int(layer_number)
                        Z_above_vec = []
                        layer_number_above_vec = []
                        Z_above_previous_vec = []
                        layer_number_above_previous_vec = []
                        theta_of_same_wetting_front_expressed_above_previous_vec = []
                        theta_of_wetting_front_below_expressed_above_previous_vec = []

                        Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3]
                        Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]
                        Z_below_previous, theta_below_previous, layer_number_below_previous = self.previous_states[l+1][0:3]
                        x_previous = psi_of_theta(theta_previous, theta_s, theta_r, n, m, alpha)
                        x_below_previous = psi_of_theta(theta_below_previous, theta_s, theta_r, n, m, alpha)

                        for k in range(0,layer_number):
                            Z_above_vec.append(max_depth_vec[k])
                            layer_number_above_vec.append(k)
                            Z_above_previous_vec.append(max_depth_vec[k])
                            layer_number_above_previous_vec.append(k)
                            theta_of_same_wetting_front_expressed_above_previous_vec.append(theta_of_psi(x_previous,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))
                            theta_of_wetting_front_below_expressed_above_previous_vec.append(theta_of_psi(x_below_previous,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))

                        prior_mass_above_vec = []
                        for k in range(0,layer_number):
                            prior_mass_above_vec.append( max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_previous_vec[k]-theta_of_wetting_front_below_expressed_above_previous_vec[k]) )
                        prior_mass_in_same_layer = Z_previous*(theta_previous-theta_below_previous)
                        prior_mass = sum(prior_mass_above_vec) + prior_mass_in_same_layer - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)

                        x = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                        x_below = psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha)
                        theta_of_same_wetting_front_expressed_above_vec = []
                        theta_of_wetting_front_below_expressed_above_vec = []
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec.append(theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))
                            theta_of_wetting_front_below_expressed_above_vec.append(theta_of_psi(x_below,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))

                        new_mass_above_vec = []
                        for k in range(0,layer_number):
                            new_mass_above_vec.append(max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k]))
                        new_mass_in_same_layer = Z*(theta-theta_below)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass>prior_mass):
                            x=x+1
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        ###15 aug
                        while (new_mass<prior_mass):
                            x=x-.1
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass>prior_mass):
                            x=x+.01
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass<prior_mass):
                            x=x-.001
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass>prior_mass):
                            x=x+.0001
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass<prior_mass):
                            x=x-.00001
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass>prior_mass):
                            x=x+.000001
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        while (new_mass<prior_mass):
                            x=x-.0000001
                            for k in range(0,layer_number):
                                theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                            theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                            for k in range(0,layer_number):
                                new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                            new_mass_in_same_layer = Z*(theta-theta_below)
                            new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                        if (l>0):
                            while (new_mass>prior_mass):
                                x=x+.00000001
                                for k in range(0,layer_number):
                                    theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                                theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                                for k in range(0,layer_number):
                                    new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-theta_of_wetting_front_below_expressed_above_vec[k])
                                new_mass_in_same_layer = Z*(theta-theta_below)
                                new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer
                        ###15 aug

                        self.current_states[l][1] = min(theta_s_vec[layer_number], theta)







                    # if (layer_number==1):
                    #     #description of mass balance and theta calculation for wetting fronts that are below the first layer
                    #     #wetting fronts in layers below the most superficial one have a more involved calculation for updating theta via mass balance. As of this writing (6 dec 2021), there is no solid documentation describing theta calculations for wetting fronts in layers below the top one. However, eventually we will produce a diagram which will demonstrate the procedure more clearly than text will. Still, a description is offered here.
                    #     #When wetting fronts are in layers below the top one, their dZdt values are fairly striaghtforward to calculate (see the derivs fucntion). However, calculation of theta via mass balance is a bit trickier. This is because each wetting front in deeper layers can be thought of as extending all the way to the surface, in terms of psi values. For example, a wetting front in layer 1 with a theta value of 0.4 will in reality extend to layer 0 with a theta value that is different (usually smaller) due to different soil hydraulic properties. But, the theta value of this extended wetting front is not recorded in states_array or current_states.
                    #     #so, simply from states_array and previous_states, the mass balance of a wetting front that, in terms of psi, extends between multiple layers cannot be calculated. Therefore, the theta values that the current wetting front *would* have in above layers is calculated from the psi value of the current wetting front, with the assumption that the hydraulic head of this wetting front is the same all the way up to the surface.
                    #     #these theta values for previous_states can be easily calculated.
                    #     #however, for theta values for the current time step, these must be iteratively calcualted via mass balance; once the new Z is calcualted for the current wetting front, a new psi value is iteratively calculated ofr the wetting front such that the new mass matches the old mass plus any relevant fluxes. Because the function that gives theta as a function of psi is highly nonlinear and is either very difficult or impossible to solve in closed form when two different soil layers are used in theta calculations, the new psi value is iteratively calculated to within a given mass balance tolerance.
                    #     #previous_states = states_array[states_array[:,4]==i]
                    #     #previous_states = states_array[states_array[:,4]==(len(states_array)-1)]
                    #     self.current_states[l][0] = self.current_states[l][0] + temp_derivs[l][0]*time_step #first the new Z value is recorded
                    #     #current_states[l][0] = current_states[l][0] + temp_derivs[l]*time_step
                    #
                    #     #then, the previous mass is calculated
                    #     Z, theta, layer_number = self.current_states[l][0:3]
                    #     layer_number = int(layer_number)
                    #     Z_above, layer_number_above = max_depth_vec[0], layer_number#-1 ###1 dec###
                    #     Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3]
                    #     Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]
                    #     Z_above_previous, layer_number_above_previous = max_depth_vec[0], layer_number-1
                    #     Z_below_previous, theta_below_previous, layer_number_below_previous = self.previous_states[l+1][0:3]
                    #     x_previous = psi_of_theta(theta_previous, theta_s, theta_r, n, m, alpha)
                    #     x_below_previous = psi_of_theta(theta_below_previous, theta_s, theta_r, n, m, alpha)
                    #     theta_of_same_wf_in_layer_above_previous = theta_of_psi(x_previous,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta_of_wf_below_in_layer_above_previous = theta_of_psi(x_below_previous,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #
                    #     prior_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above_previous-theta_of_wf_below_in_layer_above_previous)
                    #     prior_mass_in_same_layer = Z_previous*(theta_previous-theta_below_previous)
                    #     prior_mass = prior_mass_in_layer_above + prior_mass_in_same_layer - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)
                    #
                    #     x = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     x_below = psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha)
                    #     theta_of_wf_below_in_layer_above = theta_of_psi(x_below,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #
                    #     #then, the new mass is calculated. At first, this is an overestimate, because the new (increased) depth has been used with the old theta values (and relevant fluxes have been considered in prior_mass), and so new_mass is necessarily too large.
                    #     new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta-theta_below)
                    #     new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #
                    #     while (new_mass>prior_mass):
                    #         #here, increasing x makes the soil drier. This code necessarily makes x too large and therefore makes new_mass too small.
                    #         x=x+1
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         #in this next while loop however, x is nexessarily made too large again. This process repeats with successively small mass balance error tolerances until the mass balancer error becomes roughly 1e-10 mm per time step, or at least the erorr in x (psi) becomes very small and accordingly the actual mass balance error gets small too
                    #         x=x-.1
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.01
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.001
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.0001
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.00001
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.000001
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.0000001
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     # while (new_mass>prior_mass):
                    #     #     #by this point, the error in psi (x) will be very small.
                    #     #     #these successive while loops tend to be very computationally inexpensive; essentially ,we have relaced the need for dthetadt calculation, a differential equation requiring the function G (and therefore is somewhat computationally expensive) with an algebraic equation that is solved iteratively.
                    #     #     x=x+.00000001
                    #     #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     #     new_mass_in_layer_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #     #     new_mass_in_same_layer = Z*(theta-theta_below)
                    #     #     new_mass = new_mass_in_layer_above + new_mass_in_same_layer
                    #     self.current_states[l][1] = min(theta_s_vec[layer_number], theta)
                    #
                    # if (layer_number==2):
                    #     #this is the same idea for iterative calculation of the correct psi (x) value such that mass balance error is very small. Search for the text "description of mass balance and theta calculation for wetting fronts that are below the first layer" for a detailed description of how this works.
                    #     #previous_states = states_array[states_array[:,4]==i]
                    #     #previous_states = states_array[states_array[:,4]==(len(states_array)-1)]
                    #     self.current_states[l][0] = self.current_states[l][0] + temp_derivs[l][0]*time_step
                    #     #current_states[l][0] = current_states[l][0] + temp_derivs[l]*time_step
                    #
                    #     Z, theta, layer_number = self.current_states[l][0:3]
                    #     layer_number = int(layer_number)
                    #     Z_above, layer_number_above = max_depth_vec[1], layer_number#-1 ###1 dec###
                    #     Z_2_above, layer_2_number_above = max_depth_vec[0], layer_number-2
                    #     Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3]
                    #     Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]
                    #     Z_above_previous, layer_number_above_previous = max_depth_vec[1], layer_number-1
                    #     Z_2_above_previous, layer_number_2_above_previous = max_depth_vec[0], layer_number-2
                    #     Z_below_previous, theta_below_previous, layer_number_below_previous = self.previous_states[l+1][0:3]
                    #     x_previous = psi_of_theta(theta_previous, theta_s, theta_r, n, m, alpha)
                    #     x_below_previous = psi_of_theta(theta_below_previous, theta_s, theta_r, n, m, alpha)
                    #     theta_of_same_wf_in_layer_above_previous = theta_of_psi(x_previous,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta_of_wf_below_in_layer_above_previous = theta_of_psi(x_below_previous,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta_of_same_wf_in_layer_2_above_previous = theta_of_psi(x_previous,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_wf_below_in_layer_2_above_previous = theta_of_psi(x_below_previous,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #
                    #     prior_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above_previous-theta_of_wf_below_in_layer_2_above_previous)
                    #     prior_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above_previous-theta_of_wf_below_in_layer_above_previous)
                    #     prior_mass_in_same_layer = Z_previous*(theta_previous-theta_below_previous)
                    #     prior_mass = prior_mass_in_layer_2_above + prior_mass_in_layer_above + prior_mass_in_same_layer - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)
                    #
                    #     x = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     x_below = psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha)
                    #     theta_of_wf_below_in_layer_above = theta_of_psi(x_below,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta_of_wf_below_in_layer_2_above = theta_of_psi(x_below,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta-theta_below)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #
                    #     while (new_mass>prior_mass):
                    #         x=x+1
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.1
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.01
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.001
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.0001
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.00001
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass>prior_mass):
                    #         x=x+.000001
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     while (new_mass<prior_mass):
                    #         x=x-.0000001
                    #         theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #         theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #         theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #         new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #         new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #         new_mass_in_same_layer = Z*(theta-theta_below)
                    #         new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     # while (new_mass>prior_mass):
                    #     #     x=x+.00000001
                    #     #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    #     #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    #     #     new_mass_in_same_layer = Z*(theta-theta_below)
                    #     #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #     self.current_states[l][1] = min(theta_s_vec[layer_number], theta)









                if ( (l<maximum_wetting_front_number)&(layer_number_below!=layer_number) ): #in this case, we are considering wetting fronts that are not the deepest wetting front in the deepest layer but are the deepest wetting front in the layer that the wetting front occupies. Here, dZdt will be 0 from derivs, and the correct psi value has already been calculated and this information is in the theta value directly below this wetting front (new theta and Z values for wetting fronts are calculated from the bottom up, so the correct calculation has already been done), so the new theta value for this wetting front will be determined by the psi value below.
                    Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3]
                    x = psi_of_theta(theta_below, theta_s_vec[layer_number+1], theta_r_vec[layer_number+1], n_vec[layer_number+1], m_vec[layer_number+1], alpha_vec[layer_number+1])
                    theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    self.current_states[l][1] = theta

                if ( (l==maximum_wetting_front_number)&(layer_number_below!=layer_number)&(len(self.current_states)==len(parameters)) ): #this is the special case of what happens when there is only one wetting front per layer (meaning, they all have the same psi value). in this case, theta must be iteratively updated via mass balance, but dZdt = 0, free drainage boundary conditions apply so there is no capillary suction, and there are no wetting fronts below the current one.

                    #######
                    #
                    # #previous_states = states_array[states_array[:,4]==i]
                    # #previous_states = states_array[states_array[:,4]==(len(states_array)-1)]
                    #
                    # Z, theta, layer_number = self.current_states[l][0:3]
                    # layer_number = int(layer_number)
                    # Z_above, layer_number_above = max_depth_vec[1], layer_number-1
                    # Z_2_above, layer_2_number_above = max_depth_vec[0], layer_number-2
                    # Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]
                    # Z_above_previous, layer_number_above_previous = max_depth_vec[1], layer_number-1
                    # Z_2_above_previous, layer_number_2_above_previous = max_depth_vec[0], layer_number-2
                    # x_previous = psi_of_theta(theta_previous, theta_s, theta_r, n, m, alpha)
                    # theta_of_same_wf_in_layer_above_previous = theta_of_psi(x_previous,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    # theta_of_same_wf_in_layer_2_above_previous = theta_of_psi(x_previous,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #
                    # prior_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above_previous)
                    # prior_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above_previous)
                    # prior_mass_in_same_layer = Z_previous*(theta_previous)
                    # prior_mass = prior_mass_in_layer_2_above + prior_mass_in_layer_above + prior_mass_in_same_layer - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)
                    #
                    # x = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                    # theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    # theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #
                    # theta_below = 0 #these are set like this because there is no wetting front below the current one
                    # theta_of_wf_below_in_layer_above = 0
                    # theta_of_wf_below_in_layer_2_above = 0
                    #
                    # new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above-theta_of_wf_below_in_layer_2_above)
                    # new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above-theta_of_wf_below_in_layer_above)
                    # new_mass_in_same_layer = Z*(theta-theta_below)
                    # new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #
                    # while (new_mass>prior_mass):
                    #     x=x+1
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass<prior_mass):
                    #     x=x-.1
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass>prior_mass):
                    #     x=x+.01
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass<prior_mass):
                    #     x=x-.001
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass>prior_mass):
                    #     x=x+.0001
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass<prior_mass):
                    #     x=x-.00001
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass>prior_mass):
                    #     x=x+.000001
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    # while (new_mass<prior_mass):
                    #     x=x-.0000001
                    #     theta_of_same_wf_in_layer_2_above = theta_of_psi(x,theta_r_vec[layer_number-2],theta_s_vec[layer_number-2],alpha_vec[layer_number-2],n_vec[layer_number-2],m_vec[layer_number-2])
                    #     theta_of_same_wf_in_layer_above = theta_of_psi(x,theta_r_vec[layer_number-1],theta_s_vec[layer_number-1],alpha_vec[layer_number-1],n_vec[layer_number-1],m_vec[layer_number-1])
                    #     theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    #     new_mass_in_layer_2_above = max_depth_vec[0]*(theta_of_same_wf_in_layer_2_above)
                    #     new_mass_in_layer_above = max_depth_vec[1]*(theta_of_same_wf_in_layer_above)
                    #     new_mass_in_same_layer = Z*(theta)
                    #     new_mass = new_mass_in_layer_2_above + new_mass_in_layer_above + new_mass_in_same_layer
                    #
                    #
                    # self.current_states[l][1] = theta
                    # theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    # self.current_states[l][1] = min(theta_s_vec[layer_number], theta)
                    #

                    #######



                    Z, theta, layer_number = self.current_states[l][0:3]
                    layer_number = int(layer_number)
                    Z_above_vec = []
                    layer_number_above_vec = []
                    Z_above_previous_vec = []
                    layer_number_above_previous_vec = []
                    theta_of_same_wetting_front_expressed_above_previous_vec = []
                    theta_of_wetting_front_below_expressed_above_previous_vec = []

                    #Z_below, theta_below, layer_number_below = self.current_states[l+1][0:3]
                    Z_previous, theta_previous, layer_number_previous = self.previous_states[l][0:3]
                    #Z_below_previous, theta_below_previous, layer_number_below_previous = self.previous_states[l+1][0:3]
                    x_previous = psi_of_theta(theta_previous, theta_s, theta_r, n, m, alpha)
                    #x_below_previous = psi_of_theta(theta_below_previous, theta_s, theta_r, n, m, alpha)

                    for k in range(0,layer_number):
                        Z_above_vec.append(max_depth_vec[k])
                        layer_number_above_vec.append(k)
                        Z_above_previous_vec.append(max_depth_vec[k])
                        layer_number_above_previous_vec.append(k)
                        theta_of_same_wetting_front_expressed_above_previous_vec.append(theta_of_psi(x_previous,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))
                        #theta_of_wetting_front_below_expressed_above_previous_vec.append(theta_of_psi(x_below_previous,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))

                    prior_mass_above_vec = []
                    for k in range(0,layer_number):
                        prior_mass_above_vec.append( max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_previous_vec[k]-0) )
                    prior_mass_in_same_layer = Z_previous*(theta_previous-0)
                    prior_mass = sum(prior_mass_above_vec) + prior_mass_in_same_layer - (free_drainage_demand+actual_ET_demand)*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)

                    x = psi_of_theta(theta, theta_s, theta_r, n, m, alpha)
                    #x_below = psi_of_theta(theta_below, theta_s, theta_r, n, m, alpha)
                    theta_of_same_wetting_front_expressed_above_vec = []
                    #theta_of_wetting_front_below_expressed_above_vec = []
                    for k in range(0,layer_number):
                        theta_of_same_wetting_front_expressed_above_vec.append(theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))
                        #theta_of_wetting_front_below_expressed_above_vec.append(theta_of_psi(x_below,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k]))

                    new_mass_above_vec = []
                    for k in range(0,layer_number):
                        new_mass_above_vec.append(max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0))
                    new_mass_in_same_layer = Z*(theta-0)
                    new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass>prior_mass):
                        x=x+1
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    ###15 aug
                    while (new_mass<prior_mass):
                        x=x-.1
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass>prior_mass):
                        x=x+.01
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass<prior_mass):
                        x=x-.001
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass>prior_mass):
                        x=x+.0001
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass<prior_mass):
                        x=x-.00001
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass>prior_mass):
                        x=x+.000001
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer

                    while (new_mass<prior_mass):
                        x=x-.0000001
                        for k in range(0,layer_number):
                            theta_of_same_wetting_front_expressed_above_vec[k] = theta_of_psi(x,theta_r_vec[k],theta_s_vec[k],alpha_vec[k],n_vec[k],m_vec[k])
                        theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                        for k in range(0,layer_number):
                            new_mass_above_vec[k] = max_depth_vec[k]*(theta_of_same_wetting_front_expressed_above_vec[k]-0)
                        new_mass_in_same_layer = Z*(theta-0)
                        new_mass = sum(new_mass_above_vec) + new_mass_in_same_layer
                    ###15 aug

                    self.current_states[l][1] = theta
                    theta = theta_of_psi(x,theta_r_vec[layer_number],theta_s_vec[layer_number],alpha_vec[layer_number],n_vec[layer_number],m_vec[layer_number])
                    self.current_states[l][1] = min(theta_s_vec[layer_number], theta)









            #by this point, all Z and theta values have been updated. However, there are conditions that when met necessitate changing of current_states, after updating all wetting fronts via derivs and iterative psi calculation. These are described below.

            #in the event that a wetting front's theta value reaches the corresponding theta_s for its layer, that means that that wetting front has reached saturation.
            #however, the f_p calculation does not consider whether or not the soil will saturate, so the predicted infiltration would cause a theta value larger than theta_s. In the above code, if this happened, then theta was set to theta_s, which introduces a mass balance error. To correct this, the Z value of the wetting front that reached saturation is increased such that the mass balance closes (within a small tolerance). Note that the erroneously large Z value for the wetting front occurs for exactly one time step, and the amount by which Z will be too large will decrease with a smaller time step. The best way forward in the future will be to determine the exact time at which saturation would occur and split the time step at that time, but considerations involving adaptive or split time steps will happen after the FIHM, June 2022.

            #if ( (theta_previous<theta_s_vec[layer_number])&(self.current_states[l][1] == theta_s_vec[layer_number])  ):
            #if ( (self.previous_states[0][1]<theta_s_vec[0])&(self.current_states[0][1] == theta_s_vec[0])  ):
            layer_with_wf_that_supplies_free_drainage_demand = self.current_states[wf_that_supplies_free_drainage_demand][2]
            #if ( (self.previous_states[0][1]<theta_s_vec[0])&(self.current_states[wf_that_supplies_free_drainage_demand][1] == theta_s_vec[layer_with_wf_that_supplies_free_drainage_demand])  ):
            ###ok .... something I'm trying on 9 August 2022. I'm getting a small mass balance error caused I think by the first part of this if statemet, and also I don't think it's necessary. So I'm trying taking it out.
            if ( (self.current_states[wf_that_supplies_free_drainage_demand][1] == theta_s_vec[layer_with_wf_that_supplies_free_drainage_demand])  ):
            #wf_that_supplies_free_drainage_demand
                current_mass = calc_mass_bal(self.current_states)
                mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add + actual_ET_demand ### adding ET 11 feb
                #print('MBAL ERROR FOR SATURATION')
                #print(mass_balance_error)

                ###8 march test: taking out
                ##ok so indeed this fixes the mass balance error but makes the runoff look weird
                #runoff = runoff + abs(mass_balance_error)
                #actual_ET_demand = actual_ET_demand + abs(mass_balance_error)
                #
                # extra_length = mass_balance_error/(self.current_states[wf_that_supplies_free_drainage_demand][1]-self.current_states[wf_that_supplies_free_drainage_demand+1][1])
                # self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + extra_length


                #free_drainage_demand = free_drainage_demand - abs(mass_balance_error)

                ######
                while mass_balance_error<0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + .1
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand
                #
                # ###15 aug
                while mass_balance_error>0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] - .01
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error<0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + .001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error>0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] - .0001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error<0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + .00001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error>0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] - .000001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error<0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + .0000001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand

                while mass_balance_error>0:
                    self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] - .00000001
                    current_mass = calc_mass_bal(self.current_states)
                    mass_balance_error = current_mass - self.mass_vec[-1] + free_drainage_demand - precip_mass_to_add+ actual_ET_demand
                ###15 aug

                ###here's something slightly interesting: by including the next 4 lines, it is the case that the very small mass balance error introduced by overflow conditions has the opposite sign of the mass balance error intorduced by successive x (psi) calculation. This does increase the runtime in a miniscule way (like an extra 5 computations on average per time step), but it does make the mass balance error non monotonic and ultimately slightly smaller. Anyway just a small detail.
                # while mass_balance_error<0:
                #     self.current_states[wf_that_supplies_free_drainage_demand][0] = self.current_states[wf_that_supplies_free_drainage_demand][0] + .000000001
                #     current_mass = calc_mass_bal(self.current_states)
                #     mass_balance_error = current_mass - mass_vec[-1] + free_drainage_demand - precip_mass_to_add

                #######

            ###finally, this above code apparently will erroneously make the deepest wetting front too deep, in the event that it is the one satisfying the free drainage demand. In this case, the excess water will be erroneously added to free drainage, but again only for one time step and to conserve mass. The above suggested solution such that none of these small corrections are necessary involve split or adaptive time steps.
            #if ( (theta_previous<theta_s_vec[layer_number])&(current_states[l][1] == theta_s_vec[layer_number]) & (max(wetting_front_numbers)==len(parameters)) ):
            if (self.current_states[-1][0] > max_depth_vec[-1]):
                runoff_temp = runoff
                #free_drainage_demand = free_drainage_demand + (current_states[-1][0]-max_depth_vec[-1])*current_states[-1][1]#/time_step
                runoff = runoff + (self.current_states[-1][0]-max_depth_vec[-1])*self.current_states[-1][1]#/time_step
                precip_mass_to_add = precip_mass_to_add - (runoff - runoff_temp)
                self.current_states[-1][0] = max_depth_vec[-1]
                ###next try making this runoff instead



            #merging code here
            #in the event that a wetting front passes another wetting front in the same layer, a wetting front reaches the lower boundary for its layer, or in the event that a wetting front reaches the lower boundary of the deepest layer, the linked list current_states must be updated.
            for l in range(0,len(self.current_states)-0):#-0 included just because in old versions, wetting front numbering worked differently
                if (l>=max(self.wetting_front_numbers)):
                    break
                current_wf = self.current_states[l]
                current_layer = self.current_states[l][2]
                next_wf = self.current_states[l+1]
                next_layer = self.current_states[l+1][2]
                abs_depth_current = self.calc_abs_depth(current_wf)
                abs_depth_below = self.calc_abs_depth(next_wf)
                #this is if a wetting front passes the one directly below it in the same layer. In this case the wetting fronts merge.
                if ( (current_wf[0]>next_wf[0]) & (next_wf[2]==current_wf[2]) & (abs_depth_below not in boundary_depths) ):
                    current_mass_this_layer = self.current_states[l][0]*(self.current_states[l][1]-self.current_states[l+1][1]) + self.current_states[l+1][0]*(self.current_states[l+1][1]-self.current_states[l+2][1])
                    self.current_states[l][0] = current_mass_this_layer/(self.current_states[l][1]-self.current_states[l+2][1])
                    self.current_states.remove(self.current_states[l+1])
                    self.number_of_wetting_fronts = len(self.current_states)
                    #wetting_front_numbers = np.array(range(0,number_of_wetting_fronts))
                    self.wetting_front_numbers = np.arange(0,self.number_of_wetting_fronts)


            for l in range(0,len(self.current_states)-0):#-0 included just because in old versions, wetting front numbering worked differently
                if (l>=max(self.wetting_front_numbers)):
                    break
                current_wf = self.current_states[l]
                current_layer = self.current_states[l][2]
                next_wf = self.current_states[l+1]
                next_layer = self.current_states[l+1][2]
                abs_depth_current = self.calc_abs_depth(current_wf)
                abs_depth_below = self.calc_abs_depth(next_wf)
                #this is if a wetting front reaches its boundary depth for its layer.
                if ((abs_depth_current>abs_depth_below)&(abs_depth_below in boundary_depths)&(current_layer!=max_layer)): #if a wetting front reaches the boundary of its layer

                    # current_theta = min(theta_s_vec[layer_number], self.current_states[l][1])
                    # next_param_set = parameters[current_layer+1]
                    # theta_s,theta_r,n,m,alpha = theta_s_vec[current_layer],theta_r_vec[current_layer],n_vec[current_layer],m_vec[current_layer],alpha_vec[current_layer]
                    # current_psi = psi_of_theta(current_theta,theta_s,theta_r,n,m,alpha)
                    # theta_s,theta_r,n,m,alpha = theta_s_vec[current_layer+1],theta_r_vec[current_layer+1],n_vec[current_layer+1],m_vec[current_layer+1],alpha_vec[current_layer+1]
                    # new_theta = theta_of_psi(current_psi, theta_r, theta_s, alpha, n, m)
                    # overshot_mass = (self.current_states[l][1]-self.current_states[l+1][1])*(abs_depth_current-abs_depth_below)
                    # print('overshot_mass')
                    # print(overshot_mass)
                    # mbal_Z_correction = overshot_mass/current_theta
                    # self.current_states[l+1][1] = self.current_states[l][1]
                    # self.current_states.remove(self.current_states[l])
                    # self.current_states.insert(l+1,[mbal_Z_correction,min(new_theta,next_param_set[1]),current_layer+1,0,i+1]) #i+1 maybe?

                    #current_theta = min(theta_s_vec[layer_number], self.current_states[l][1])
                    current_theta = min(theta_s_vec[current_layer], self.current_states[l][1])
                    overshot_depth = abs_depth_current-abs_depth_below
                    old_bdy_theta = self.current_states[l+1][1]
                    current_param_set = parameters[current_layer]
                    next_param_set = parameters[current_layer+1]
                    current_theta_s, current_theta_r, current_n, current_m, current_alpha = current_param_set[1], current_param_set[0], current_param_set[4], current_param_set[5], current_param_set[3]
                    next_theta_s, next_theta_r, next_n, next_m, next_alpha = next_param_set[1], next_param_set[0], next_param_set[4], next_param_set[5], next_param_set[3]
                    current_psi_value = psi_of_theta(current_theta, current_theta_s, current_theta_r, current_n, current_m, current_alpha)
                    new_theta_value = theta_of_psi(current_psi_value,next_theta_r,next_theta_s,next_alpha,next_n,next_m)
                    mbal_correction = overshot_depth*(current_theta-old_bdy_theta)
                    mbal_Z_correction = mbal_correction/(new_theta_value-self.current_states[l+2][1])
                    self.current_states[l+1][1] = self.current_states[l][1]
                    self.current_states.remove(self.current_states[l])
                    self.current_states.insert(l+1,[mbal_Z_correction,min(new_theta_value,next_param_set[1]),current_layer+1,0,i])


                    #current_states.insert(l+1,[0.000000001*0+mbal_Z_correction,min(new_theta_value+0*0.000000001,next_param_set[1]),current_layer+1]) #leaving this line in as a comment, perhaps superstitiously, just in the event that I eventually get mysterious errors related to division by zero once a wetting front passes the lower boundary of its layer



            #in the event that a wetting front passes another wetting front in the same layer, a wetting front reaches the lower boundary for its layer, or in the event that a wetting front reaches the lower boundary of the deepest layer, the linked list current_states must be updated.
            for l in range(0,len(self.current_states)-0):#-0 included just because in old versions, wetting front numbering worked differently
                if (l>=max(self.wetting_front_numbers)):
                    break
                current_wf = self.current_states[l]
                current_layer = self.current_states[l][2]
                next_wf = self.current_states[l+1]
                next_layer = self.current_states[l+1][2]
                abs_depth_current = self.calc_abs_depth(current_wf)
                abs_depth_below = self.calc_abs_depth(next_wf)
                #this is if a wetting front passes the one directly below it in the same layer. In this case the wetting fronts merge.
                if ( (current_wf[0]>next_wf[0]) & (next_wf[2]==current_wf[2]) & (abs_depth_below not in boundary_depths) ):
                    current_mass_this_layer = self.current_states[l][0]*(self.current_states[l][1]-self.current_states[l+1][1]) + self.current_states[l+1][0]*(self.current_states[l+1][1]-self.current_states[l+2][1])
                    self.current_states[l][0] = current_mass_this_layer/(self.current_states[l][1]-self.current_states[l+2][1])
                    self.current_states.remove(self.current_states[l+1])
                    self.number_of_wetting_fronts = len(self.current_states)
                    #wetting_front_numbers = np.array(range(0,number_of_wetting_fronts))
                    self.wetting_front_numbers = np.arange(0,self.number_of_wetting_fronts)


                #this is what happens when a wetting front reaches the bottom of the deepest layer.
                #in this case, the free drainage is made erroneosuly large for one time step. This is because technically, when the maximum depth for the soil profile is achieved, it is actually slightly exceeded. The small amount of water that exceeded the maximum depth is added to the free drainage, so the free drainage is erroneously high for a time step but mass balancer closes. This is similar to the recently described solution to infiltration that would cause theta>theta_s; the amount by which free drainage is made erroneously large is made smaller when the time step is smaller, and a good solution would be to split the time step such that the maximum depth is exactly met.

            for l in range(0,len(self.current_states)-0):#-0 included just because in old versions, wetting front numbering worked differently
                if (l>=max(self.wetting_front_numbers)):
                    break
                current_wf = self.current_states[l]
                current_layer = self.current_states[l][2]
                next_wf = self.current_states[l+1]
                next_layer = self.current_states[l+1][2]
                abs_depth_current = self.calc_abs_depth(current_wf)
                abs_depth_below = self.calc_abs_depth(next_wf)

                if ((abs_depth_current>abs_depth_below)&(abs_depth_below in boundary_depths)&(current_layer==max_layer)):

                    #free_drainage_demand = free_drainage_demand + (current_states[l][0]-current_states[l+1][0])*(current_states[l][1]-current_states[l+1][1])
                    ###23 dec replaced current_states[l+1][0] with max_depth_vec[-1] due to an error where sometimes the deepest wetting front will exceed its max possible depth slightly
                    free_drainage_demand = free_drainage_demand + (self.current_states[l][0]-max_depth_vec[-1])*(self.current_states[l][1]-self.current_states[l+1][1])
                    #free_drainage_augmented = 1

                    current_theta = self.current_states[l][1]
                    self.current_states[l+1][1] = current_theta
                    self.current_states.remove(self.current_states[l])
                    total_num_wetting_fronts = len(self.current_states)
                    #wetting_front_numbers = np.array(range(0,total_num_wetting_fronts))
                    self.wetting_front_numbers = np.arange(0,total_num_wetting_fronts)

                    #3 feb takeout
                    #number_of_wetting_fronts = len(current_states)
                    #wetting_front_numbers = np.array(range(0,number_of_wetting_fronts))

            ### here, 23 dec 2021, and 27 dec, trying for merging code in the event that the top WF becomes too dry via free drainage or ET

            for l in range(0,len(self.current_states)-0):#-0 included just because in old versions, wetting front numbering worked differently
                if (l>=max(self.wetting_front_numbers)):
                    break
                current_wf = self.current_states[l]
                current_layer = self.current_states[l][2]
                next_wf = self.current_states[l+1]
                next_layer = self.current_states[l+1][2]
                abs_depth_current = self.calc_abs_depth(current_wf)
                abs_depth_below = self.calc_abs_depth(next_wf)


                if ( (current_wf[1]<(next_wf[1])) & (current_layer==next_layer) ): #
                    #print('free drainage (or ET) was large enough to make a wetting front erroneously have less moisture than the WF below it in the same layer.')
                    # print('current_wf')
                    # print(current_wf)
                    # print('next_wf')
                    # print(next_wf)

                    self.current_states.remove(self.current_states[l])
                    self.number_of_wetting_fronts = len(self.current_states)
                    self.wetting_front_numbers = np.arange(0,self.number_of_wetting_fronts)

                    #ok. so previous approach was to correct the actual ET itself.
                    #the approach I think I should do now is to remove the most superficial wetting front and then distribute the original actual ET demand throughout the new wf that contributes to ET, which I calculate




                    # print(' ')
                    # print('MBAL CORRECTED FOR ET')
                    # print(' ')

                    # if (current_layer==0):
                    #     theta = self.current_states[l][1]
                    #     self.current_states[0][1] = theta



                    if (current_layer > 0):
                        theta = self.current_states[l][1]
                        theta_s = theta_s_vec[current_layer]
                        theta_r = theta_r_vec[current_layer]
                        n = n_vec[current_layer]
                        m = m_vec[current_layer]
                        alpha = alpha_vec[current_layer]
                        x = psi_of_theta(theta,theta_s,theta_r,n,m,alpha)

                        for k in range(0,current_layer):
                            theta_s = theta_s_vec[k]
                            theta_r = theta_r_vec[k]
                            n = n_vec[k]
                            m = m_vec[k]
                            alpha = alpha_vec[k]
                            self.current_states[k][1] = theta_of_psi(x, theta_r, theta_s, alpha, n, m)

                        # theta_s = theta_s_vec[1]
                        # theta_r = theta_r_vec[1]
                        # n = n_vec[1]
                        # m = m_vec[1]
                        # alpha = alpha_vec[1]
                        # self.current_states[1][1] = theta_of_psi(x, theta_r, theta_s, alpha, n, m)


                    # if (current_layer==1):
                    #     theta = self.current_states[l][1]
                    #     theta_s = theta_s_vec[current_layer]
                    #     theta_r = theta_r_vec[current_layer]
                    #     n = n_vec[current_layer]
                    #     m = m_vec[current_layer]
                    #     alpha = alpha_vec[current_layer]
                    #     x = psi_of_theta(theta,theta_s,theta_r,n,m,alpha)
                    #
                    #     theta_s = theta_s_vec[0]
                    #     theta_r = theta_r_vec[0]
                    #     n = n_vec[0]
                    #     m = m_vec[0]
                    #     alpha = alpha_vec[0]
                    #     self.current_states[0][1] = theta_of_psi(x, theta_r, theta_s, alpha, n, m)
                    #
                    # if (current_layer==2):
                    #     theta = self.current_states[l][1]
                    #     theta_s = theta_s_vec[current_layer]
                    #     theta_r = theta_r_vec[current_layer]
                    #     n = n_vec[current_layer]
                    #     m = m_vec[current_layer]
                    #     alpha = alpha_vec[current_layer]
                    #     x = psi_of_theta(theta,theta_s,theta_r,n,m,alpha)
                    #
                    #     theta_s = theta_s_vec[0]
                    #     theta_r = theta_r_vec[0]
                    #     n = n_vec[0]
                    #     m = m_vec[0]
                    #     alpha = alpha_vec[0]
                    #     self.current_states[0][1] = theta_of_psi(x, theta_r, theta_s, alpha, n, m)
                    #
                    #     theta_s = theta_s_vec[1]
                    #     theta_r = theta_r_vec[1]
                    #     n = n_vec[1]
                    #     m = m_vec[1]
                    #     alpha = alpha_vec[1]
                    #     self.current_states[1][1] = theta_of_psi(x, theta_r, theta_s, alpha, n, m)




                    ### 5 July takeout and move to after creation of new wetting front
                    ###ET_caused_merge = True
                    previous_h_p = self.h_p_vec[-1]
                    change_in_mass_temp = calc_mass_bal(self.current_states) - calc_mass_bal(self.previous_states) + self.h_p - previous_h_p
                    #actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff)
                    ###actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff)
                    if (i>0):
                        if ((precip_data[i]>0)&(precip_data[i-1]==0)):
                            actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_mass_to_add-0*precip_data[i]*time_step+runoff)
                        else:
                            actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_mass_to_add-precip_data[i]*time_step+runoff)
                    ###

                    # if (change_in_mass_temp>0):
                    #     actual_ET_demand = precip_to_add_for_mbal

                    #21 march, trying *0 with precip_to_add_for_mbal

                    #mbal_error_temp = change_in_mass_temp + free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff+actual_ET_demand
                    #actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-precip_to_add_for_mbal-0*precip_data[i]*time_step+0*runoff)

            #here a new wetting front is created if there is precip in the current time step but no precip in the previous time step
            #also note that a new wetting front won't be created at time = 0 if there is nonzero precip at time = 0. This is due to the fact that we don't know the precip data before the first time step so I make the assumption that precip at time = 0 cintributes to an already existing wetting front.
            if ( (i>0) & (self.h_p==0) ):
                #here, a new wetting front will be created if it's not the first time step, and if there was no precip in the previous time step but there is precip in the current time step.
                #the initial depth of this new wetting front is calculated using the concept of "dry depth", described in the 2015 GARTO paper (An efficient and guaranteed stable numerical method for continuous modeling of infiltration and redistribution with a shallow dynamic water table).
                if ((precip_data[i]>0)&(precip_data[i-1]==0)&(self.current_states[0][1]>=theta_s_vec[0])):
                    runoff = precip_data[i]*time_step #+runoff*0 #in the event that we have completely saturated soil with 0 ET, which could happen for example when using a PET method that relies on solar radiation at night with multiple intense precip pulses
                if ((precip_data[i]>0)&(precip_data[i-1]==0)&(self.current_states[0][1]<theta_s_vec[0])):
                    temp_theta = self.current_states[0][1]
                    tau = time_step*K_s_vec[0]/(theta_s_vec[0]-temp_theta)
                    psi_b = self.psi_b_vec[0]
                    lambdaa = self.lambda_vec[0]
                    if self.closed_form_capillary_drive_term:
                        G_new_wf = G_closed(theta_s_vec[0],temp_theta,psi_b,lambdaa,theta_r_vec[0],theta_s_vec[0])
                    else:
                        G_new_wf = G(psi_of_theta(theta_s_vec[0], theta_s_vec[0], theta_r_vec[0], n_vec[0], m_vec[0], alpha_vec[0]), psi_of_theta(temp_theta, theta_s_vec[0], theta_r_vec[0], n_vec[0],m_vec[0],alpha_vec[0]), alpha_vec[0], n_vec[0], m_vec[0], K_s_vec[0])
                    dry_depth = min( 0.5*(tau + ( (tau**2) + 4*tau*G_new_wf )**0.5 ) , max_depth_vec[0] ) ### coefficient "x" in x*(tau should be 0.5. For basic sensitivity analysis to dry depth, sometimes I change it, but the theory from 2015 GAR paper says it's 0.5
                    ###!!! note that dry depth originally has a factor of 0.5 in front
                    if (dry_depth==max_depth_vec[0]):
                        print('dry depth greater than layer 0 thickness, so dry depth has been set to layer 0 thickness')
                        #I think it is unlikely that this will ever happen; still, if it does, IDK if the mass balancer will work out as currently implemented
                    new_theta = min(temp_theta + precip_data[i]*time_step/dry_depth, theta_s_vec[0])
                    if (new_theta==theta_s_vec[0]):
                        #hm interesting, the dry depth concept assumes that, for the first time step of a wetting front's existence, runoff is only possible as overflow runoff. That is, infiltration capacity excess runoff is not possible. I suppose this makes sense because f_p calculations require G which requires a difference in theta between the current wetting front and the wetting front below in the same layer. At the very start of precip, this difference is 0, and it's in the denominator, so f_p can't be calculated for the very first time step ... this will probably only lead to small errors, especially with a small time step. And frankly, with an adaptive time step, the problem will be as small as you want it to be via error tolerances.
                        print('initial precip was intense enough to fill up dry depth. So gonna need to simulate runoff here.')
                        runoff = runoff*0 + precip_data[i]*time_step - (theta_s_vec[0]-temp_theta)*dry_depth #March 18: "runoff +"" part necessary?
                        precip_mass_to_add = precip_data[i]*time_step - runoff
                        #self.h_p = self.h_p + precip_data[i]*time_step - (theta_s_vec[0]-temp_theta)*dry_depth

                    precip_mass_to_add = precip_data[i]*time_step - runoff + self.h_p_vec[-1]
                    self.current_states.insert(0,[dry_depth,new_theta,0,0,i])
                    self.number_of_wetting_fronts = len(self.current_states)
                    #wetting_front_numbers = np.array(range(0,number_of_wetting_fronts))
                    self.wetting_front_numbers = np.arange(0,self.number_of_wetting_fronts)

            # previous_h_p = self.h_p_vec[-1]
            # change_in_mass_temp = calc_mass_bal(self.current_states) - calc_mass_bal(self.previous_states) + self.h_p - previous_h_p
            # #actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff)
            # actual_ET_demand = -1*(change_in_mass_temp + free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff)


            #so if runoff has been predicted, and ponded depth is less than its max value, the runoff will go to ponding rather than running off
            if (runoff>0):
                if (self.h_p==params.h_p_max):
                    runoff=runoff

                else:
                    if ((runoff+self.h_p)>self.h_p_max):
                        h_p_temp_run = self.h_p
                        self.h_p = self.h_p_max
                        #runoff = runoff - (self.h_p_max - self.h_p_vec[-1])
                        runoff = runoff - (self.h_p_max - h_p_temp_run)
                    else:
                        if (ok_to_increment_h_p):
                            self.h_p = self.h_p + runoff
                            runoff = 0
                        else:
                            self.h_p = self.h_p


            for wf_number_set in range(0,len(self.current_states)): #this just sets wetting front number and time step to be the right values in the linked list representing the current states
                self.current_states[wf_number_set][3]=wf_number_set
                self.current_states[wf_number_set][4]=i



            self.current_states[-1][0] = self.max_depth_vec[-1]



            ### 29 dec code that fixes error where sometimes at end of precip event that causes saturation, actual infiltration is erroneously high ###
            # if ( (current_states[0][1]==theta_s_vec[0]) & ( calc_mass_bal(current_states)+precip_mass_to_add - (runoff + free_drainage_demand) > np.dot(theta_s_vec,max_depth_vec) ) ):
            #     print('the thing')
            #     diff = (calc_mass_bal(current_states)+precip_mass_to_add - (runoff + free_drainage_demand )) - np.dot(theta_s_vec,max_depth_vec)
            #     precip_mass_to_add = precip_mass_to_add - diff
            #     runoff = runoff + diff


            #previous_states = current_states

            self.actual_ET_vec.append(actual_ET_demand)

            #free_drainage_vec.append(free_drainage_demand)

            self.h_p_vec.append(self.h_p)
            self.runoff_vec.append(runoff/time_step)

            #runoff_vec[i] = runoff/time_step

            self.actual_infiltration_vec.append(precip_mass_to_add/time_step)
            #actual_infiltration_vec[i] = precip_mass_to_add/time_step

            #finally, states_array, which records past values of current_states, is updated, along with wetting front numbers and a time step.
            if i in time_steps_to_record_profile:
                current_states_temp = np.array(self.current_states)
                #current_states_temp = np.column_stack([current_states_temp, wetting_front_numbers])
                #current_states_temp = np.column_stack([current_states_temp, i+1+np.zeros(number_of_wetting_fronts)])
                self.states_array = np.vstack([self.states_array,(current_states_temp) ])

            self.wetting_front_depths = np.array(self.current_states)[:,0]
            self.wetting_front_moistures = np.array(self.current_states)[:,1]

            self.precip_mm_per_h = precip_data[i]
            # print('the twilight zone')
            # print(self.wetting_front_moistures)
            self.PET_mm_per_h = PET_data[i]
            self.runoff_mm_per_h = runoff/time_step
            self.actual_infiltration_mm_per_h = precip_mass_to_add/time_step
            self.bottom_flux_mm_per_h = free_drainage_demand/time_step
            self.actual_ET_mm_per_h = actual_ET_demand/time_step
            self.ponding_depth_mm = self.h_p

            # print('actual ET demand at end of time step (mm/h)')
            # print(actual_ET_demand/time_step)
            # print(' ')

            #just used for mass balance monitoring
            precip_to_add_for_mbal = precip_mass_to_add
            self.fluxes.append(free_drainage_demand-0*precip_to_add_for_mbal-precip_data[i]*time_step+runoff+actual_ET_demand)

            previous_h_p = self.h_p_vec[-2]
            change_in_h_p = previous_h_p - self.h_p


            self.h_p_fluxes.append(change_in_h_p)
            #fluxes_cumulative.append(sum(fluxes))

            self.free_drainage_flux_vec.append(free_drainage_demand/time_step)
            #free_drainage_flux_vec[i] = free_drainage_demand/time_step


            change_in_mass = calc_mass_bal(self.current_states)-calc_mass_bal(self.previous_states)-change_in_h_p

            self.mass_bal_vec.append(change_in_mass+self.fluxes[-1]+0*self.h_p_fluxes[-1])
            self.mass_balance_error_value_cumulative = self.mass_balance_error_value_cumulative + change_in_mass+self.fluxes[-1]+0*self.h_p_fluxes[-1]

            #else:
            #    mass_bal_vec.append(mass_bal_vec[-1]+change_in_mass+fluxes[-1])
            #print(len(mass_bal_vec))


            ### these are for BMI/CSDMS conventions
            # self.precip_mm_per_h = precip_data[i]
            # self.PET_mm_per_h = PET_data[i]
            # self.runoff_mm_per_h = runoff/time_step
            # self.actual_infiltration_mm_per_h = precip_mass_to_add/time_step
            # self.bottom_flux_mm_per_h = free_drainage_demand/time_step
            # self.actual_ET_mm_per_h = actual_ET_demand/time_step
            # self.ponding_depth_mm = self.h_p
            #self.wetting_front_depths #already defined
            #self.wetting_front_moistures #already defined

            if self.first_time_step == True:
                self.fluxes = [self.fluxes.pop(0)]
                self.runoff_vec = [self.runoff_vec.pop(0)]
                self.actual_infiltration_vec = [self.actual_infiltration_vec.pop(0)]
                #self.f_p=0
                self.actual_ET_vec = [self.actual_ET_vec.pop(0)]
                self.mass_vec = np.delete(self.mass_vec, 0)
                self.free_drainage_flux_vec = [self.free_drainage_flux_vec.pop(0)]
                self.mass_bal_vec = [self.mass_bal_vec.pop(0)]
                #self.h_p_vec = [self.h_p_vec.pop(0)]

                    # print(' ')
                    # print('TEST FOR ERROR IN NUMBER OF WETTING FRONTS')
                    # print(' ')

            self.first_time_step = False

            if self.verbose:
                print(' ')
                print(' ')
                print('time step:')
                print(i)
                print('precipitation (mm/h):')
                print(self.precip_data[i])
                print('Linked list representing all wetting fronts at this time step:')
                print(self.current_states)

            # if self.verbose:
            #     print('theta at which AET=0.5*PET')
            #     print(theta_50)
            #     print('corresponding psi value')
            #     print(h_50)
            if self.verbose:
                print(' ')
                print(' ')
                print('time step:')
                print(i)
                print('precipitation (mm/h):')
                print(self.precip_data[i])
                #print(self.current_states)
                print('current number of wetting fronts:')
                print(max(self.wetting_front_numbers)+1)
                print('previous mass in soil (mm):')
                print(calc_mass_bal(self.previous_states))
                print('previous h_p value:')
                print(self.h_p_vec[-1])
                print('previous states:')
                print(self.previous_states)
                print('current mass in soil (mm):')
                print(calc_mass_bal(self.current_states))
                print('h_p (mm):')
                print(self.h_p)
                print('current mass in system (h_p + soil water):')
                print(calc_mass_bal(self.current_states)+self.h_p)
                print('current states:')
                print(self.current_states)
                print('mass difference (new-old, mm):')
                print(calc_mass_bal(self.current_states)-calc_mass_bal(self.previous_states)+self.h_p-previous_h_p)
                print('net mass of water leaving soil (mm):')
                print(self.fluxes[-1])
                print('h_p difference (new-old, mm):')
                print(self.h_p-previous_h_p)
                print('MASS BALANCE ERROR for this time step (mm):')
                print(self.mass_bal_vec[-1])
                print('cumulative mass balance error (mm):')
                print(self.mass_balance_error_value_cumulative)
                print('actual ET:')
                print(actual_ET_demand)
                print('bottom flux:')
                print(free_drainage_demand)
                print('infiltration to soil:')
                print(precip_mass_to_add)
                print('mass added to system (including h_p):')
                print(precip_data[i]*time_step-runoff)
                print('runoff:')
                print(runoff)
                print('f_p:')
                print(f_p)
                print('psi values for wetting fronts at end of time step (-mm):')
                for WF in self.current_states:
                    theta_s,theta_r,n,m,alpha = theta_s_vec[WF[2]],theta_r_vec[WF[2]],n_vec[WF[2]],m_vec[WF[2]],alpha_vec[WF[2]]
                    psi_val_to_print = psi_of_theta(WF[1],theta_s,theta_r,n,m,alpha)
                    print(psi_val_to_print)
                print(' ')
                print(' ')

            self.current_time = time_step*(i+1)
            self.current_time_step = self.current_time_step+1






    # #this is just because a different name is used for plotting later, or used to be anyway
    # actual_infil_vec = actual_infiltration_vec









    # def finalize(self):
    #     self.runtime_1 = datetime.now() #this records the end time for computations
    #     np.save('states_array',self.states_array)
    #     print("run success! computation runtime:")
    #     print(self.runtime_1-self.runtime_0)


    # forcing_data = forcing_data[0:len(precip_data)]
    # forcing_data['runoff[mm/h]'] = runoff_vec
    # forcing_data['actual_infil[mm/h]'] = actual_infil_vec
    # h_p_vec = h_p_vec[:len(precip_data)]
    # forcing_data['ponded_head[mm]'] = h_p_vec
    # forcing_data['bottom_flux[mm/h]'] = free_drainage_flux_vec
    # forcing_data['bottom_demand[mm]'] = free_drainage_vec
    # forcing_data['water_in_soil[mm]'] = mass_vec
    # forcing_data['cumulative_mbal_error(mm)'] = mass_balance_error_vec
    # forcing_data['actual_ET_per_step(mm)'] = actual_ET_vec
    #
    # forcing_data.to_pickle("LGAR_output.pkl")






# import config as config
# length_of_simulation = config.length_of_simulation
#
# test=LGAR('config_2')
# test.run_model(length_of_simulation)
# #test.finalize()
# print('this text means that the class ran successfully')
#
# # print('before update')
# # print(test.current_states)
# #test.update_until(3)
# # test.update()
# # print('after update')
# # print(test.current_states)
#
# test.initalize('config_2')
# print('init test')
# print(test.current_states)
