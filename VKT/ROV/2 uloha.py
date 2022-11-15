# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import lmfit
from inspect import signature
from matplotlib.ticker import AutoMinorLocator

path='C:/Users/aleks/Projects/arameic-mishmash/VKT/olej_rot/'

odtek=pd.read_csv(f'{path}CmCvm.csv', header=0, delimiter=';', decimal=',')

def plot_param(ax): #not mine
    FONT = 'Times New Roman'   
    plt.ylabel('C[ml/s]', fontname = FONT, fontweight = 'bold', fontsize = 18)
    plt.xlabel('p [Pa]',  fontname = FONT, fontweight = 'bold', fontsize = 18)
    # plt.xlim(0.,0.5)   #!!!!!
    # plt.ylim(0.,0.5)
    #graphic parameters
    plt.xticks(fontname=FONT, fontweight = 'bold', fontsize = 15)
    plt.yticks(fontname=FONT, fontweight = 'bold', fontsize = 15)
    
    ax.tick_params(which = 'major', direction='out', length=6, width=1.5)
    ax.tick_params(which = 'minor', direction='out', length=3, width=1)
    
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(1.5)
        
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    
    ax.grid(which = 'major', c = 'gray', linewidth = 0.5, linestyle = 'solid') 
    ax.grid(which = 'minor', c = 'gray', linewidth = 0.3, linestyle = 'dashed')             
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    leg = plt.legend(  loc = 'best', shadow = True, fancybox=False)
    
    leg.get_frame().set_linewidth(1)
    leg.get_frame().set_edgecolor('k')
    
    for text in leg.get_texts():
          plt.setp(text, fontname=FONT, fontsize = 14)


def print_params(x, y, func, p0=1, maxfev=1000, bounds=(-np.inf, np.inf)):
    #fit with curve_fit and print parameters
    par, cov = curve_fit(func, x, y, p0=p0,  maxfev = maxfev, bounds=bounds)
    standev=np.sqrt(np.diag(cov))
    for i, val in enumerate(par):
        print(f'{i+1}. parameter = {val:.2e} ± {standev[i]:.2e}')
    print('%%%%%%')
    return(par, standev)


def split_fit(data, fitfunction, N=1, maxfev=1000,\
              param_hint=None, split_index=None, fit_smooth=1000, bounds=None, extraplus=0, extraminus=0):
    #split data weather to N pieces, or by split_index and fit each separately 
    #space - skip number of points between fits
    #param_hint 2-d array with initial tips for params example:
    #param_hint = np.array([[1,2,3],[3,4,5]]) for 2 subfits with three parameters
    #split_index - 2-d array of row indexes where to split and fit
    #bounds example: 2 subint, 3 pars: [[[xmin1,xmin2,xmin3],
    #                                    [xmax1,xmax2,xmax3]]
    #                                   [[ymin1,ymin2,ymin3],
    #                                    [ymax1,ymax2,ymax3]]]
    #extraplus/minus - will just draw fit lines outside data point (units are same as x axes)
    if type(split_index)!=type(None):
        split_index.sort() #in case if they're not in ascending order
    split_data=pd.DataFrame()
    #split the data we use in fit
    if split_index==None:
        spl = np.array_split(data, N)         #split the data we want to fit
        for i in range(N):                    #for cycle is needed becouse of syntax used later
            split_data=pd.concat((split_data, pd.concat({i:spl[i]})))
    else:                     
        for i, subindex in enumerate(split_index):
            split_data = pd.concat((split_data, pd.concat({i :\
                                    data.iloc[split_index[i][0]:split_index[i][1]]})))
    if type(bounds)==type(None): #np.shape(bounds)=(#subintervals, 2, #pars)
        bounds=np.full((split_data.index[-1][0]+1,2,\
                        len(signature(fitfunc).parameters)-1), np.inf)
        bounds[:,0]*=-1
    subfits = pd.DataFrame()
    if type(param_hint)==type(None):
        param_hint=np.full(split_data.index[-1][0]+1, None)
    for sub in range(split_data.index[-1][0]+1): #finds number of subsets and for each subset makes a fit
        pars, standevs = print_params(split_data.loc[sub].iloc[:, 0],\
                                      split_data.loc[sub].iloc[:, 1],\
                                      fitfunction, p0=param_hint[sub],\
                                      maxfev = maxfev, bounds=bounds[sub])
        xaxes = np.linspace(split_data.loc[sub].iloc[0, 0]-extraminus, \
                          split_data.loc[sub].iloc[-1, 0]+extraplus, fit_smooth)
        yaxes = fitfunction(xaxes, *pars)
        subfits = pd.concat((subfits, pd.concat({sub :pd.DataFrame([xaxes,yaxes]).T})))        
        a = subfits.copy() # !!! Despair IDK why can't I just do subfits.columns=data.columns...
        a.columns=data.columns
    return(a)
    
def fitfunc(x, a, b):
    y=a*x+b
    return y

cejch=odtek

cejch_fit0=split_fit(cejch[['p_str','C']][3:],fitfunc, extraplus=10, extraminus=10)
cejch_fit1=split_fit(cejch[['p_str','Cvm']],fitfunc, extraplus=10, extraminus=25)

fig, ax=plt.subplots()

dataline=ax.scatter(cejch['p_str'], cejch['C'] ,color='red', s=15, marker="D")
dataline2=ax.scatter(cejch['p_str'], cejch['Cvm'] ,color='green', s=15, marker="P")


fit=ax.plot(cejch_fit0['p_str'], cejch_fit0['C'] ,color='red',\
            label=r'$C$',alpha=0.7, linewidth=2)
fit1=ax.plot(cejch_fit1['p_str'], cejch_fit1['Cvm'] ,color='green',\
            label=r'$C_{vm}$',alpha=0.7, linewidth=2)

            
plot_param(ax)
    