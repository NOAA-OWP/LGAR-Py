###This .py file reformats raw USDA SCAN data to the format used by LGAR-Py.
###only the following 4 lines should be changed by the user
forcing_data_folder = '~/desktop/LGAR-py/forcing_data_files/Fort_Assiniboine/' #This is the location of both the raw and reformatted forcing datasets.
raw_forcing_data_file_name = 'forcing_data_Fort_Assiniboine_raw.csv' #This is the input (raw) USDA SCAN forcing dataset.
formatted_forcing_data_name = 'forcing_data_resampled_Fort_Assiniboine.csv' #this will be the name of the output forcing data, in the correct format for LGAR-Py.
time_step_formatting = 300/3600 #this is the output time step of the resampled forcing data, expressed in hours. The default value of 300/3600 is 5 minutes expressed in hours.






def forcing_data_formatter_fxn(path_string,raw_forcing_data_file_name,formatted_forcing_data_name,freq):

    import pandas as pd
    import numpy as np

    raw_forcing_data_file_string = path_string + raw_forcing_data_file_name
    forcing_data = pd.read_csv(raw_forcing_data_file_string,header=0, delimiter=',')

    alpha = 1.3
    rho_w = 1000

    def e_s(T):
        return(611*np.exp((17.27*T)/(273.3+T)))

    def delta(T):
        return((4098*e_s(T))/((273.3+T)**2))

    def l_v(T):
        return(2500-2.36*T)

    def gamma(T):
        return(66.8)
        #return(1.005*101325*pressure/(0.622*l_v(T)))
        #note that gamma is often given as the constant 66.8 pascals/decree C. the other version, commented out, takes temperature and pressure (input pressure in atm, and temp as deg C) as inputs.

    def E_r(R_n,T):
        return(R_n/(l_v(T)*1000*rho_w))

    def E_PT(R_n,T):
        return(E_r(R_n,T)*alpha*delta(T)/(delta(T)+gamma(T)))

    def F_to_C(F):
        return( (F - 32) * 5.0/9.0 )


    PET_vec = []
    precip_vec = []

    for i in (range(0,len(forcing_data))):
        PET_vec.append(3.6e6*E_PT(float(forcing_data['Solar Radiation Average (watt/m2)'][i]),F_to_C(float(forcing_data['Air Temperature Average (degF)'][i]))))
        precip_vec.append(float(forcing_data['Precipitation Increment (in)'][i])*25.4)

    forcing_data['PET(mm/h)'] = PET_vec
    forcing_data['P(mm/h)'] = precip_vec

    forcing_data = forcing_data.filter(['P(mm/h)', 'PET(mm/h)','Date'])

    df = pd.to_datetime(forcing_data['Date'])

    forcing_data['Datetime_date'] = df

    forcing_data.set_index('Datetime_date', inplace=True)

    resampling_frequency = str(int(freq*3600))+'S'
    forcing_data = forcing_data.resample(resampling_frequency).interpolate()


    forcing_data = forcing_data.filter(['P(mm/h)', 'PET(mm/h)'])
    #
    forcing_data.to_csv(path_string + formatted_forcing_data_name) #this is for saving forcing data / sending to HYDRUS

    return(forcing_data)


forcing_data_formatter_fxn(path_string=forcing_data_folder, raw_forcing_data_file_name=raw_forcing_data_file_name, formatted_forcing_data_name=formatted_forcing_data_name,freq=time_step_formatting)
