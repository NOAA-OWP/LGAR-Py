def forcing_data_formatter_fxn(path_string,raw_forcing_data_file_name,formatted_forcing_data_name,freq):


    import pandas as pd
    import numpy as np

    raw_forcing_data_file_string = path_string + raw_forcing_data_file_name

    forcing_data = pd.read_csv(raw_forcing_data_file_string,header=0, delimiter=',')
    #forcing_data = pd.read_csv("forcing_data/Phillipsburg_data/Phillipsburg_data.csv",header=0, delimiter=',')

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
        PET_vec.append(3.6e6*E_PT(float(forcing_data['Phillipsburg (2093) Solar Radiation Average (watt/m2)'][i]),F_to_C(float(forcing_data['Phillipsburg (2093) Air Temperature Average (degF)'][i]))))
        precip_vec.append(float(forcing_data['Phillipsburg (2093) Precipitation Increment (in)'][i])*25.4)

    forcing_data['PET(mm/h)'] = PET_vec
    forcing_data['P(mm/h)'] = precip_vec

    forcing_data = forcing_data.filter(['P(mm/h)', 'PET(mm/h)','Date'])

    df = pd.to_datetime(forcing_data['Date'])

    forcing_data['Datetime_date'] = df

    forcing_data.set_index('Datetime_date', inplace=True)

    resampling_frequency = str(int(freq*3600))+'S'
    forcing_data = forcing_data.resample(resampling_frequency).interpolate()

    ##########
    #this creates synthetic forcing data for testing or concept visualization purposes
    #for i in range(0,len(forcing_data)):

    #for i in range(0,90000):
    for i in range(0,1000):
        forcing_data['P(mm/h)'][i] = forcing_data['P(mm/h)'][i]*0
        forcing_data['PET(mm/h)'][i] = forcing_data['PET(mm/h)'][i]*0

    precip_number = 7#40

    for i in range(0,6):
        forcing_data['P(mm/h)'][i] = 0
    for i in range(6,11):
        forcing_data['P(mm/h)'][i] = precip_number
    for i in range(11,59):
        forcing_data['P(mm/h)'][i] = 0*precip_number
    for i in range(59,99+30):
        forcing_data['P(mm/h)'][i] = precip_number
    for i in range(99+30,120):
        forcing_data['P(mm/h)'][i] = 0*precip_number

    # for i in range(400,410):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(410,420):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(420,430):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(430,440):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(440,450):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(450,460):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(460,470):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(470,480):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(480,490):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(490,500):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(500,510):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(510,520):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(520,530):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(530,540):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(540,550):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(550,560):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(560,570):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(570,580):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(580,590):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(590,600):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(600,610):
    #     forcing_data['P(mm/h)'][i] = 0

    # for i in range(0,200):
    #     forcing_data['P(mm/h)'][i]=forcing_data['P(mm/h)'][i+400]
    #
    # for i in range(400,600):
    #     forcing_data['P(mm/h)'][i]=forcing_data['P(mm/h)'][i+400]





    #
    # for i in range(0,50):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(10,20):
    #     forcing_data['P(mm/h)'][i] = 10
    # for i in range(50,100):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(100,150):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(150,200):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(200,250):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(250,300):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(300,350):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(350,400):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(400,450):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(450,500):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(500,550):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(550,600):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(600,650):
    #     forcing_data['P(mm/h)'][i] = precip_number*0
    # for i in range(650,700):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(700,750):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(750,800):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(800,850):
    #     forcing_data['P(mm/h)'][i] = precip_number*3
    # for i in range(850,900):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(900,950):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(950,1000):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1000,1050):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1050,1100):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(1100,1150):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1150,1200):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1200,1250):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1250,1300):
    #     forcing_data['P(mm/h)'][i] = 0
    # for i in range(1300,1350):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1350,1400):
    #     forcing_data['P(mm/h)'][i] = precip_number
    # for i in range(1400,1450):
    #     forcing_data['P(mm/h)'][i] = precip_number


    # def sinusoidal_synthetic_PET(step):
    #     result = -0.4*np.cos(2*np.pi/288*step)
    #     if result<0:
    #         result = 0
    #     return(result)
    #
    # for i in range(0,1500):
    #     forcing_data['PET(mm/h)'][i] = sinusoidal_synthetic_PET(i)

    #########






    forcing_data = forcing_data.filter(['P(mm/h)', 'PET(mm/h)'])
    return(forcing_data)
    #
    #forcing_data.to_csv('forcing_data_resampled.csv') #this is for saving forcing data / sending to HYDRUS
