# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:15:28 2019

@author: ?
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from matplotlib import cm

def assemble(para, cache):
    """ Assemble linear system Jacobian * dx = F
    
    Process:
        0. Obtain relevant informations
        1. Loop over grid:
            1.1 Deal with BC node at x=0
            1.2 Deal with BC node at x=L
            1.3 Deal with interior nodes
            1.4 Obtain values on imaginary nodes (Ug1 and Ug2)
                for two BCs
            1.4 Assemble Jacobian (a diagonal matrix)
        2. Calculate temperature gradient dT2
        3. Assemble F
    
    Return: dictionary containing cache data
    """


    # BC informations
    typeX0 = para['x=0 type']
    valueX0 = para['x=0 value']
    typeXL = para['x=L type']
    valueXL = para['x=L value']
    ne = cache['ne']
    Z = cache['Zbar']
    Kn = cache['Kn']
    Kb = para['boltzman']
    x = para['x']
    dx = para['deltaX']

    numberOfNode = int(para['numberOfNode'])
    # Containers
    T = cache['T']; T0 = cache['T0']        #let T=T[i,j] then T0=T[i, j-1]
    F = cache['F']; Jacobian = cache['Jacobian']
    dt = cache['dt']  


    alphas, betas, heatflux = cache['alpha'], cache['beta'], cache['heatflux']
        ##Coulomb logarithm 
    coulog = 23-np.log(np.sqrt(ne)*Z/T**1.5) #np.ones(len(para['x'])) #23-np.log(np.sqrt(ne)*Z/T**1.5)
    
        ##Thermal velocity (profile)
    v=np.sqrt(T*Kb/para['m_e'])
        ##Lambda mean free path
    lamb = v**4/(ne*para['Gamma']*coulog)*1/np.sqrt(Z+1)
    gradT=np.gradient(T,x)

    ##Knudsen number according to (5)
    Kn = -lamb*gradT/T
    kappa = para['conductivity']*1.31e10/coulog*para['tau']**(cache['beta']-5/2)
    cache['kappa_LOCAL'] = para['conductivity']*1.31e10/coulog*para['tau']
    heatflux = -(para['SpitzerHarmCOND']/Z)*((Z+0.24)/(Z+4.2))*T**2.5*gradT




    '''    
    Loop over grid
    '''


    for i in range(0, numberOfNode):
        # BC node at x=0
        if i == 0:
            if typeX0 == 'heatFlux':
                Ug1 = fixedGradient(valueX0, kappa[i], dx, T[0], alphas[i], betas[i]) #boundary values
                Jacobian[0][1] = (1/dx**2)*(alphas[i+1]*kappa[i+1]*betas[i+1]*T[i+1]**(betas[i+1]-1)*T[i]\
                                           -(betas[i+1]+1)*alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] - alphas[i]*kappa[i]*T[i]**betas[i])
            elif typeX0 == 'fixedTemperature':
                Ug1 = fixedValue(valueX0, T[1])
                Jacobian[0][1] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] + alphas[i]*kappa[i]*Ug1**betas[i]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(Ug1+T[i+1]))
        # BC node at x=L
        elif i == numberOfNode-1:
            if typeXL == 'heatFlux':
                Ug2 = fixedGradient(valueXL, kappa[i], dx, T[-1], alphas[i], betas[i])  #boundary values
                Jacobian[-1][-2] = (1/dx**2)*(betas[i-1]*kappa[i-1]*alphas[i-1]*T[i-1]**(betas[i-1]-1)*T[i]\
                                           -(betas[i-1]+1)*alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1] - alphas[i]*kappa[i]*T[i]**betas[i])
            elif typeXL == 'fixedTemperature':
                Ug2 = fixedValue(valueXL, T[-2])
                Jacobian[-1][-2] = 0
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i]*kappa[i]*Ug2**betas[i] + alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(T[i-1]+Ug2))  
        # Interior nodes

        else:   #!!! \alpha_{i+1/2} := alpha[i]
            Jacobian[i][i+1] = (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*betas[i+1]*T[i+1]**(betas[i+1]-1)*T[i]\
                                           -(betas[i+1]+1)*alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] - alphas[i]*kappa[i]*T[i]**betas[i])
            
            Jacobian[i][i-1] = (1/dx**2)*.5*(betas[i-1]*kappa[i-1]*alphas[i-1]*T[i-1]**(betas[i-1]-1)*T[i]\
                                           -(betas[i-1]+1)*alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1] - alphas[i]*kappa[i]*T[i]**betas[i])
            
            Jacobian[i][i] = (3/2*ne[i]*Kb)/dt + (1/dx**2)*.5*(alphas[i+1]*kappa[i+1]*T[i+1]**betas[i+1] + alphas[i-1]*kappa[i-1]*T[i-1]**betas[i-1]\
                            +2*(betas[i]+1)*alphas[i]*kappa[i]*T[i]**betas[i] - (betas[i]*alphas[i]*kappa[i]*T[i]**(betas[i]-1))*(T[i-1]+T[i+1]))


    # Calculate F (right hand side vector)
    d2T = secondOrder(T, Ug1, Ug2, alphas, betas,kappa)
    F = (3/2*ne)*(T - T0)*Kb/dt + d2T/dx**2 # Vectorization   dT/dt - a d2T/dx2=F/dt

    # Store in cache
    cache['F'] = F; cache['Jacobian'] = Jacobian
    cache['kappa'], cache['Kn'], cache['heatflux'] = kappa, Kn, heatflux
    cache['coulog'] = coulog
    return cache


def initialize(para):
    """ Initialize key data
    
    T: current step temperature
    T0: last step temperature
    TProfile: temperature results in time and space
    F: B as right hand side of Ax = B
    Jacobian: A as left had side of Ax = B
    
    Return: a dictionary
    """
    
    numberOfNode = int(para['numberOfNode'])
    numOfTimeStep = para['numberOfTimeStep']
    T_init = para['InitTeProfile']
    Zbar_init=para['InitZbarProfile']
    ne_init=para['InitneProfile']
    Kn_init=para['InitKnProfile']
    heatflux_init= para['InitHeatflux']
    T = T_init 
    T0 = T_init
    alpha_init = para['alphas']
    beta_init = para['betas']
    


    #Define empty matrices that will contain time evolution of the profiles
    TProfile = np.zeros((numberOfNode, numOfTimeStep + 1))
    heatflux_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Zbar_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    Kn_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    ne_prof = np.zeros((numberOfNode, numOfTimeStep + 1))
    F = np.zeros((numberOfNode, 1))
    Jacobian = np.zeros((numberOfNode, numberOfNode))

    #Filling first column with initial values of the quantities
    TProfile[:,0] = T.reshape(1,-1)
    heatflux_prof[:,0] = heatflux_init.reshape(1,-1)
    Kn_prof[:,0] = Kn_init.reshape(1,-1)
    ne_prof[:,0] = ne_init.reshape(1,-1)
    Zbar_prof[:,0] = Zbar_init.reshape(1,-1)

    times=np.array([0])

    coulog_init = 23-np.log(np.sqrt(ne_init)*Zbar_init/T_init**1.5)

    dt = Exception("dt wasn't calculated")
    kappa = Exception("kappa wasn't calculated")
    cache = {'T':T,'T0':T0,'TProfile':TProfile, 'alpha':alpha_init, 'beta':beta_init, 'heatflux':heatflux_init,
             'F':F,'Jacobian':Jacobian, 'time':0, 'times':times, 'dt':dt, 'kappa': kappa, 'Zbar':Zbar_init, 
             'ne':ne_init,'Kn':Kn_init, 'Kn_prof':Kn_prof,'ne_prof':ne_prof,'Zbar_prof':Zbar_prof,
             'Log':pd.DataFrame(),'heatflux_prof':heatflux_prof, 'coulog':coulog_init}
    return cache


def solveLinearSystem(para, cache):
    """ Solve Ax=B
    
    Process:
        1. Get A = Jacobian matrix (Jacobian)
        2. Get B = Right hand side equation (F)
        3. Calculate dT
        4. Update T
        5. Store in cache
        
    Return: a dictionary
    """
    relax = para['relaxation'] 
    A = cache['Jacobian']
    B = cache['F']
    dT = np.linalg.solve(A, B)
    T = cache['T']
    T = T-dT * relax       #T(j+1)=T(j)+JI``(F)
    T[np.where(T<=0)]=10
    cache['T'] = T
    cache['dT'] = dT
    return cache

def storeUpdateResult(cache):
    """ Store results
    Update T0
    Store temperaure results into a dataframe and 
    save it in the cache.
    """
    
    timeStep = cache['ts']
    TProfile = cache['TProfile'] 
    heatflux_prof = cache['heatflux_prof']    
    Zbar_prof = cache['Zbar_prof']   
    Kn_prof = cache['Kn_prof']   
    ne_prof = cache['ne_prof']


    alpha = cache['alpha']      #current profile
    beta = cache['beta']
    heatflux = cache['heatflux']
    Zbar = cache['Zbar']
    Kn = cache['Kn']
    ne = cache['ne']
    T = cache['T']
    cache['T0'] = T.copy()

    TProfile[:,timeStep] = T.reshape(1,-1)
    heatflux_prof[:,timeStep] = heatflux.reshape(1,-1)
    Zbar_prof[:,timeStep] = Zbar.reshape(1,-1)
    Kn_prof[:,timeStep] = Kn.reshape(1,-1)
    ne_prof[:,timeStep] = ne.reshape(1,-1)

    return cache

def newtonIteration(para, cache):
    """ Newton's Iteration for Equation System
    
    Process:
        1. Get max iteration, convergence limit
        2. Call assemble function to get Jacobian and F(RHS)
        3. Solve for dT, update solution
        4. Evaluate F, get value of 2-norm
        5. If solution converged, break, output to screen and
           return cache.
    
    """
    
    maxIteration = para['maxIteration']
    convergence = para['convergence']

    T = cache['T'];     #let T=T[i,j] then T0=T[i, j-1] 
    #cache['dt'] = para['Time_multiplier']*np.min(3/2*para['InitneProfile']*para['boltzman']*para['deltaX']**2/((para['conductivity']*1.31e10/cache['coulog']*para['tau']**(cache['beta']-5/2))*cache['alpha']*T**2.5))
    cache['dt'] = para['dt']
    cache['time']+=cache['dt']
    cache['times'] = np.append(cache['times'],cache['time'])
    log = cache['Log']
    ts = cache['ts']




    if para['Break_condition']=='max_iter':
        for n in range(maxIteration):
            cache = assemble(para, cache)
            F = cache['F']
            norm = np.linalg.norm(F)
            energy=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
            
            if n==0: slump, energy_init = np.copy(norm), np.copy(energy)
            #if norm/np.linalg.norm(cache['T']) < convergence:
            elif np.abs(energy_init-energy)/energy_init < convergence and n!=0:
                log.loc[ts,'PhysicalTime'] = cache['time']
                log.loc[ts,'Iteration'] = n+1
                log.loc[ts,'Residual'] = np.abs(energy_init-energy)/energy_init
                break
            cache = solveLinearSystem(para, cache)

    elif para['Break_condition']=='lower_bound':
        n=0
        while np.abs(energy_init-energy)/energy_init >= convergence:
            cache = assemble(para, cache)
            F = cache['F']
            norm = np.linalg.norm(F)
            energy=np.mean(1.5*para['boltzman']*cache['T']*cache['ne'])
            if n==0: slump, energy_init = np.copy(norm), np.copy(energy)
            cache = solveLinearSystem(para, cache)
            n+=1
    else: 
        print('Wrong break condition')
        quit()

    print('[{:3.0f}'.format(ts), ']',
          '[{:6.2E}'.format(cache['time']),']',
          '[{:2.0f}'.format(n+1), ']',
          #'[{:8.2E}'.format(norm/energy),']',
          '[{:8.2E}'.format(np.abs(energy_init-energy)/energy_init),']', #shows what portion of energy was lost in one time step
          '[{:8.2E}'.format(norm/slump),']',
          '[{:8.2E}'.format(np.max(cache['beta'])),']',
          '[{:8.2E}'.format(np.max(cache['alpha'])),']',
          '[{:8.2E}'.format(np.min(T)),']',
          '[{:8.2E}'.format(np.max(T)),']',
          #' [','{:8.2E}'.format(np.mean(cache['T'])),']')
          '[{:8.16E}'.format(energy),']')
    cache['Log']=log #update Log
    return cache


def solve(para, FluxLimiter=None):
    """ Main function to solve heat conduction
    
    Input: a Pandas series containing all parameters
    
    Process:
        1. Initialize cache
        2. Time marching 
        3. Newton's iteration for discretized PDE for singe time 
           step
        4. Update T, save result to T profile
    
    Return: temperature profile as final result
    """
    
    print(" Heat Conduction Solver")
    start = time.time()
    cache = initialize(para)
    cache['FluxLimiter']=FluxLimiter
    numOfTimeStep = para['numberOfTimeStep']
    print(' [Step] [Time] [Iter] [Residue] [Newton outcome] [Max beta] [Max alpha] [Minimal T] [Maximal T] [meanEnergy]')
    for timeStep in range(1, numOfTimeStep+1):
        cache['ts'] = timeStep
        cache = newtonIteration(para, cache)
        cache = storeUpdateResult(cache)
    TProfile = pd.DataFrame(cache['TProfile'], columns=cache['times'],index=para['x'])
    heatflux_prof = cache['heatflux_prof']
    runtime = time.time() - start
    print('[Cost] CPU time spent','%.3f'%runtime,'s')
    return TProfile, cache, heatflux_prof


"Boundary condition"

def fixedValue(value, U2):
    """  Dirichlet boundary condition

    Assume that value of variable at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    
    Ug = 2 * value - U2 
    return Ug


def fixedGradient(q, kappa, dx, U1, alphas, betas):
    """  Neumann boundary condition
    
    Assume that the resulted gradient at BC is fixed.
    Please see any numerical analysis text book for details.
    
    Return: float
    """
    #(U1-Ug)/dx*kappa=q
    Ug =  q *(betas+1)/(kappa*alphas) * 2 * dx  + U1
    return Ug



def secondOrder(U, Ug1, Ug2, alphas, betas, kappa):
    """ Calculate second order derivative
    
    Centered differencing approximation.
    D2U/Dx2 = (U[i-1] - 2U[i] + U[i+1])/dx**2
    
    For BC nodes, use the values on ghost nodes.
    
    Ug1: value on ghost node at x=0
    Ug2: value on ghost node at x=L
    
    Please see any numerical analysis text book for details.
    
    Return: numpy array
    """
    
    d2U = np.zeros((U.size,))

    for i in range(0, U.size):
        if i==0:
            d2U[i] = .5*(alphas[i]*kappa[i]*Ug1**betas[i] + alphas[i]*kappa[i]*U[i]**betas[i])*(U[i] - Ug1)\
                    -.5*(alphas[i+1]*kappa[i+1]*U[i+1]**betas[i+1] + alphas[i]*kappa[i]*U[i]**betas[i])*(U[i+1] - U[i])
        elif i==(U.size - 1):
            d2U[i] = .5*(alphas[i-1]*kappa[i-1]*U[i-1]**betas[i-1] + alphas[i]*kappa[i]*U[i]**betas[i])*(U[i] - U[i-1])\
                    -.5*(alphas[i]*kappa[i]*Ug2**betas[i] + alphas[i]*kappa[i]*U[i]**betas[i])*(Ug2 - U[i])
        else:
            d2U[i] = .5*(alphas[i-1]*kappa[i-1]*U[i-1]**betas[i-1] + alphas[i]*kappa[i]*U[i]**betas[i])*(U[i] - U[i-1])\
                    -.5*(alphas[i+1]*kappa[i+1]*U[i+1]**betas[i+1] + alphas[i]*kappa[i]*U[i]**betas[i])*(U[i+1] - U[i])

    return d2U

#visualisation
def evolutionField(results, name):
    """ Generate 3D temperature fields
    
    For better understanding of the results
    
    Inputs:
        1. parameter, a pandas series
        2. results, a numpy array
    """
    
    X = results.index
    Y = results.columns*1e9
    X, Y = np.meshgrid(X, Y)
    
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel(r'$x$ [cm]', fontsize=16,labelpad=15)
    ax.set_ylabel(r'$t$ [ns]', fontsize=16,labelpad=15)
    ax.set_zlabel(name, fontsize=16,labelpad=15)


    ax.grid(visible=None, which='minor', axis='both')
    Z = results.T.values
    ax.plot_surface(X, Y, Z, 
                    cmap=cm.seismic,
                    linewidth=0, 
                    antialiased=True)
    plt.show()
