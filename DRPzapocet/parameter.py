# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 22:41:22 2022

@author: ?
"""
import heatConduction as hc
import pandas as pd
import numpy as np


def main():
    """ Generate parameter
    
    1. Generate system-level parameters
    2. Generate material properties, grid, time, bcs
    
    Return: a Pandas series
    """
    
    #### Import initial profile (x, Te, ne, Zbar)
    column = 'values'
    df = pd.Series(name=column, dtype='float64')
    df = df.astype('object')
    
    # Grid
    length = 1000
    df.at['numberOfNode'] = length
    df.at['x'] = np.linspace(0, 1, length)
    df.at['length'] = df['x'][-1]
    
    # Initial conditions
    df.at['InitTeProfile'] = (np.sin((df['x']-df['x'].min())/(df['length']-df['x'].min())*6*np.pi)**2 + 1.4)
    df.at['InitneProfile'] = np.linspace(1, 1,length)
    df.at['alphas'] =  np.linspace(1, 1,length) 
    df.at['betas'] = np.linspace(2.5, 0, length)
    
    # Material
    df.at['material function'] = 'Given by ne, alpha and beta'
    df.at['boltzman'] = 1 #1.602178e-12 #8.617333262e-5 (eV/K)   1.38e-16 (CGS)

    # Boundary conditions
    df.at['x=0 type'] = 'heatFlux' #'heatFlux' - Neumann or 'fixedTemperature' - Dirichlet
    df.at['x=0 value'] = 0
    df.at['x=L type'] = 'heatFlux' #'heatFlux' - Neumann or 'fixedTemperature' - Dirichlet
    df.at['x=L value'] = 0
    
   
    # Solution
    df.at['Time_multiplier'] = 100
    df.at['deltaX'] = df['x'][11] - df['x'][10]
    df.at['dt'] = df['Time_multiplier']*np.mean(3/2*df['InitneProfile']*df['boltzman']*df['deltaX']**2/\
                               (df['alphas']*df['InitTeProfile']**2.5))
    df.at['Break_condition'] = 'max_iter' #'max_iter'/'lower_bound'   #Chooses what condition will stop newton iteration 
    df.at['MaxTime'] = 200
    df.at['numberOfTimeStep'] = 100 #int(df['MaxTime']/df['dt'])
    df.at['maxIteration'] = 30
    df.at['convergence'] = 1e-9
    df.at['relaxation'] = 1 #value in [0-1] Very sensitive!!!


    df.at['InitgradTeProfile'] = np.gradient(df['InitTeProfile'], df['x'])
    df.at['InitHeatflux'] = -df['alphas'] * df['InitTeProfile']**2.5*df['InitgradTeProfile']
    return df



if __name__ == "__main__":
    parameter = main()
    results, cache = hc.solve(parameter)
    T = pd.DataFrame(results)

    hc.evolutionField(T[T.columns[:]], r'$T$')
    
    
    