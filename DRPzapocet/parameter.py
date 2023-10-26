# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: ?
"""
import heatConduction as hc
import pandas as pd
import numpy as np
import scipy.constants as const





def main():
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    #### Import initial profile (x, Te, ne, Zbar)
    init_profile=pd.read_csv('./DRPzapocet/initial_profile.csv', index_col=(0))#/DRPzapocet


    # init_profile=pd.DataFrame(np.transpose(np.array([init_profile['x'], np.linspace(2000, 1000, 400), \
    #                                                  np.linspace(init_profile['ne'][176],init_profile['ne'][176],400), \
    #                                                 np.linspace(init_profile['Zbar'][176],init_profile['Zbar'][176],400)])),
    #                           columns=init_profile.columns,index=init_profile.index)


    column = 'values'
    df = pd.Series(name = column, dtype='float64')
    df = df.astype('object')
    
    # Grid
    df.at['Time_multiplier'] = 1e-8
    df.at['length'] = init_profile['x'].iloc[-1]
    df.at['numberOfNode'] = int(len(init_profile))
    df.at['x']=init_profile['x'].values
    
    #hrani si
    length=1000 #len(init_profile)
    df.at['numberOfNode'] = length
    df.at['x']=np.linspace(init_profile['x'].min(), init_profile['x'].max(), length)
    teplota=1000*(np.sin((df['x']-df['x'].min())/(df['length']-df['x'].min())*6*np.pi)**2+1.4)
    # Initial conditions
    df.at['InitTeProfile'] = teplota#init_profile['Te'].values
    df.at['InitneProfile'] = np.linspace(init_profile['ne'].mean(),init_profile['ne'].mean(), length)#init_profile['ne'].values
    df.at['InitZbarProfile'] = np.linspace(init_profile['Zbar'].mean(),init_profile['Zbar'].mean(), length)#init_profile['Zbar'].values

    alpha=np.linspace(1,1,length)
    alpha[198:202]=1
    beta=np.linspace(2.5,2.5, length)
    beta[200:]=0
    df.at['alphas']= alpha#np.linspace(1,1, len(init_profile))#np.linspace(1,8, len(init_profile))
    df.at['betas'] = np.linspace(2.5,0, length)#np.linspace(0,2.5, len(init_profile)) 
    
    # Material
    df.at['material function'] = 'Given by ne, Z, alpha and beta'
    df.at['conductivity'] = (df['InitZbarProfile']+0.24)/(df['InitZbarProfile']*(df['InitZbarProfile']+4.2))
    df.at['tau'] = 1e-3 #look for eq (4) in  Calculation of Heat Conduction Utilizing Neural Networks
    df.at['boltzman']=1.602178e-12#8.617333262e-5 (eV/K)   1.38e-16 (CGS)
    df.at['m_e'] = 9.1094*1e-28 #g (CGS)
    df.at['q_e'] = 4.8032*1e-10 #cm3/2 g1/2 s-1 (CGS)
    df.at['Gamma'] = 4 * const.pi * df['q_e']**4/df['m_e']**2
    df.at['SpitzerHarmCOND'] = 6.1e+02

    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0
    
   
    # Solution
    df.at['deltaX'] = df['x'][11]-df['x'][10]  #for different [i] dx differs at 16th decimal place
    df.at['dt']=df['Time_multiplier']*np.mean(3/2*df['InitneProfile']*df['boltzman']*df['deltaX']**2/\
                               (df['conductivity']*df['alphas']*df['InitTeProfile']**2.5))
    df.at['Break_condition'] = 'max_iter' #'max_iter'/'lower_bound'   #Chooses what condition will stop newton iteration 
    df.at['MaxTime'] = 2e-11
    df.at['numberOfTimeStep'] =int(df['MaxTime']/df['dt'])

    df.at['maxIteration'] = 30
    df.at['convergence'] = 1e-9
    df.at['relaxation'] =1# value in [0-1] Very sensitive!!!

    df.at['InitgradTeProfile'] = np.gradient(df['InitTeProfile'], df['x'])

    ##Knudsen number according to (5)
    lamb = np.sqrt(df['InitTeProfile']*df['boltzman']/df['m_e'])**4/\
        (df['InitneProfile']*df['Gamma']*23-np.log(np.sqrt(df['InitneProfile'])\
        *df['InitZbarProfile']/df['InitTeProfile']**1.5))*1/np.sqrt(df['InitZbarProfile']+1)
    
    df.at['InitKnProfile'] = -lamb*df['InitgradTeProfile']/df['InitTeProfile']
    df.at['InitHeatflux'] = -(df['SpitzerHarmCOND']/df['InitZbarProfile'])*((df['InitZbarProfile']+0.24)/\
                            (df['InitZbarProfile']+4.2))*df['InitTeProfile']**2.5*df['InitgradTeProfile']
    return df



if __name__ == "__main__":
    

    parameter = main()
    results, cache, heatflux = hc.solve(parameter)
    T = pd.DataFrame(results)

    heatflux3d=pd.DataFrame(heatflux, columns=T.columns,index=T.index)
    hc.evolutionField(T[T.columns[:]], r'$T$ [eV]')
    
    
    