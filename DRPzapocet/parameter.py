# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: ?
"""
import postprocessing as pp
import heatConduction as hc
import pandas as pd
import numpy as np
import scipy.constants as const
#### Import initial profile, used in PyTorch part. (x, Te, gradTe, ne, Zbar)
init_profile=pd.read_csv('./Data/init_profile.csv', index_col=(0))#/DRPzapocet
####
init_profile=init_profile.iloc[::100,:]
init_profile.reset_index(drop=True, inplace=True)
#init_profile['Te']/=1.001**((init_profile['Te']))
#init_profile['Te']=1000
#init_profile['Te'].iloc[100:180]=1000
#init_profile['ne']/=init_profile['ne']*3/2
#####


def main():
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    column = 'values'
    df = pd.Series(name = column, dtype='float64')
    df = df.astype('object')
    
    # Grid
    df.at['Time_multiplier'] = 3e-7
    df.at['length'] = init_profile['x'].iloc[-1]
    df.at['numberOfNode'] = int(len(init_profile))
    df.at['x']=init_profile['x'].values

    # Material
    df.at['material function'] = 'Given by ne, Z, alpha and beta'
    df.at['conductivity'] = (init_profile['Zbar']+0.24)/(init_profile['Zbar']*(init_profile['Zbar']+4.2)).values
    df.at['tau'] = 1e-3 #look for eq (4) in  Calculation of Heat Conduction Utilizing Neural Networks
    df.at['boltzman']=1.602178e-12#8.617333262e-5 (eV/K)   1.38e-16 (CGS)
    df.at['m_e'] = 9.1094*1e-28 #g (CGS)
    df.at['q_e'] = 4.8032*1e-10 #cm3/2 g1/2 s-1 (CGS)
    df.at['Gamma'] = 4 * const.pi * df['q_e']**4/df['m_e']**2

    # Initial conditions
    df.at['InitTeProfile'] = init_profile['Te'].values
    df.at['InitneProfile'] = init_profile['ne'].values
    df.at['InitgradTeProfile'] = init_profile['gradTe'].values
    df.at['InitZbarProfile'] = init_profile['Zbar'].values
    df.at['InitKnProfile'] = init_profile['Kn'].values

    
    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0
    
    #NN

    df.at['alphas']= np.linspace(1,1, len(init_profile)) #precal_alpha #np.linspace(1,8, len(init_profile))
    df.at['betas'] = np.linspace(2.5,2.5, len(init_profile)) #precal_beta #np.linspace(2.5,2.5, len(init_profile)) 
    df.at['heatflux'] = np.linspace(0,0, len(init_profile))

   
    # Solution
    df.at['deltaX'] = df['x'][11]-df['x'][10]  #for different [i] dx differs at 16th decimal place
    df.at['dt']=df['Time_multiplier']*np.min(3/2*df['InitneProfile']*df['boltzman']*df['deltaX']**2/\
                               (df['conductivity']*df['alphas']*df['InitTeProfile']**2.5))
    df.at['Break_condition'] = 'max_iter' #'max_iter'/'lower_bound'   #Chooses what condition will stop newton iteration 
    df.at['MaxTime'] = 4e-10
    df.at['numberOfTimeStep'] = int(df['MaxTime']/df['dt']) #Is replaced by MaxTime

    df.at['maxIteration'] = 30
    df.at['convergence'] = 1e-9
    df.at['relaxation'] =1# value in [0-1] Very sensitive!!!




    return df



if __name__ == "__main__":
    

    parameter = main()
    results, cache, alphas, betas, heatflux = hc.solve(parameter)
    T = pp.preprocess(parameter, results)
    pp.evolutionField(T)
    #np.linspace(parameter['x'][0], parameter['x'].iloc[-1], 8 )   #0-L  TODO: global variable?
    positions = T.index[::int(len(init_profile['x'])*3e-2)]
    pp.thermalCouplePlot(T, positions)
    times = T.columns[::int(len(T.columns)/10)][1:4]
        #'numberOfTimeStep'*'deltaTime'  TODO: global variable?
    pp.temperatureDistribution(T, times)
    
    
    
    