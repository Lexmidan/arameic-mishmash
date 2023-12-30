# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2019

@author: ?
"""
import heatConduction as hc
import pandas as pd
import numpy as np
import os
from pathlib import Path


def main(path):
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    #### Import initial profile (x, Te, ne, Zbar)
    init_profile = pd.read_csv(f'{path}/initial_profile.csv', index_col=(0))
    column = 'values'
    df = pd.Series(name=column, dtype='float64')
    df = df.astype('object')
    
    # Grid
    df.at['Time_multiplier'] = 100
    df.at['length'] = init_profile['x'].iloc[-1]
    df.at['numberOfNode'] = int(len(init_profile))
    df.at['x'] = init_profile['x'].values
    
    #hrani si
    length=1000 #len(init_profile)
    df.at['numberOfNode'] = length
    df.at['x'] = np.linspace(init_profile['x'].min(), init_profile['x'].max(), length)


    teplota=(np.sin((df['x']-df['x'].min())/(df['length']-df['x'].min())*6*np.pi)**2+1.4)
    # Initial conditions
    df.at['InitTeProfile'] = teplota#init_profile['Te'].values
    df.at['InitneProfile'] = np.linspace(1, 1,length)#np.linspace(init_profile['ne'].mean(),init_profile['ne'].mean(), length)#init_profile['ne'].values

    alpha = np.linspace(1, 1,length)

    beta = np.linspace(2.5, 2.5, length)
    beta[200:]=0

    df.at['alphas'] = alpha  #np.linspace(1,1, len(init_profile))#np.linspace(1,8, len(init_profile))
    df.at['betas'] = np.linspace(2.5, 0, length) #np.linspace(0,2.5, len(init_profile)) 
    
    # Material
    df.at['material function'] = 'Given by ne, alpha and beta'
    df.at['boltzman'] = 1#.602178e-12 #8.617333262e-5 (eV/K)   1.38e-16 (CGS)

    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' or 'fixedTemperature'
    df.at['x=L value'] = 0
    
   
    # Solution
    df.at['deltaX'] = df['x'][11] - df['x'][10]  #for different [i] dx differs at 16th decimal place
    df.at['dt'] = df['Time_multiplier']*np.mean(3/2*df['InitneProfile']*df['boltzman']*df['deltaX']**2/\
                               (df['alphas']*df['InitTeProfile']**2.5))
    df.at['Break_condition'] = 'max_iter' #'max_iter'/'lower_bound'   #Chooses what condition will stop newton iteration 
    df.at['MaxTime'] = 200
    df.at['numberOfTimeStep'] = 100#int(df['MaxTime']/df['dt'])

    df.at['maxIteration'] = 30
    df.at['convergence'] = 1e-9
    df.at['relaxation'] = 1 # value in [0-1] Very sensitive!!!

    df.at['InitgradTeProfile'] = np.gradient(df['InitTeProfile'], df['x'])

    ##Knudsen number according to (5)
    df.at['InitHeatflux'] = -df['alphas'] * df['InitTeProfile']**2.5*df['InitgradTeProfile']
    return df



if __name__ == "__main__":
    path = Path(os.getcwd())

    parameter = main(path)
    results, cache, heatflux = hc.solve(parameter)
    T = pd.DataFrame(results)

    heatflux3d = pd.DataFrame(heatflux, columns=T.columns, index=T.index)
    hc.evolutionField(T[T.columns[:]], r'$T$ [eV]')
    
    
    