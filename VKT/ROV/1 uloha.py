import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import lmfit
from inspect import signature

path='C:/Users/aleks/Projects/arameic-mishmash/VKT/olej_rot/'

bez_prop=pd.read_csv(f"{path}bez_prop.csv", header=0)
bez_prop['p1[Pa]']=bez_prop['p1[Pa]'].str.replace(',','.').astype(float)
s_prop=pd.read_csv(f'{path}1_s_prop - List 1.csv', header=0)


def fitfunc(x, SV):
    p=np.exp(np.log(96e3)-SV*(x))
    return p
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
fit_s_prop=split_fit(s_prop[:9], fitfunc, maxfev=10000,\
                         param_hint=np.array([[1e-2]]), \
                         bounds=np.array([[[1e-4],[3e-2]]]),\
                         extraplus=13)
fit_bez_prop=split_fit(bez_prop[:11], fitfunc, maxfev=10000,\
                         param_hint=np.array([[1e-2]]), \
                         bounds=np.array([[[1e-4],[3e-2]]]),\
                         extraplus=18)
# x=np.linspace(0,900,1000)
# y=fitfunc(x, par[0])c
fig1, axs1 = plt.subplots(2, 2, figsize=(5, 5))
axs1[0][0].semilogy(s_prop['t[s]'], s_prop['p2[Pa]'], label="S proplachovanim")
axs1[0][0].semilogy(fit_s_prop['t[s]'][0], fit_s_prop['p2[Pa]'][0],\
                 label="Fit s prop.1", ls='--',alpha=0.8, color='red')


axs1[1][0].semilogy(bez_prop['t[s]'], bez_prop['p1[Pa]'], label="Bez proplachovani")
axs1[1][0].semilogy(fit_bez_prop['t[s]'][0],fit_bez_prop['p1[Pa]'][0],\
                 label="Fit bez prop.1", ls='--', alpha=0.8, color='red')
# axs1[1].semilogy(fit_bez_prop['t[s]'][1],fit_bez_prop['p1[Pa]'][1],\
#                  label="Fit bez prop.2", ls='--',alpha=0.8, color='red')
axs1[0][1].plot(s_prop['t[s]'], s_prop['p2[Pa]'], label="S proplachovanim")
axs1[0][1].plot(fit_s_prop['t[s]'][0], fit_s_prop['p2[Pa]'][0],\
                 label="Fit s prop.1", ls='--',alpha=0.8, color='red')
axs1[1][1].plot(bez_prop['t[s]'], bez_prop['p1[Pa]'], label="Bez proplachovani")
axs1[1][1].plot(fit_bez_prop['t[s]'][0],fit_bez_prop['p1[Pa]'][0],\
                 label="Fit bez prop.1", ls='--', alpha=0.8, color='red')

axs1[1][1].set_xlabel(r'Čas $[\mathrm{s}]$')
axs1[1][0].set_xlabel(r'Čas $[\mathrm{s}]$')
axs1[1][0].set_ylabel(r'Tlak $[\mathrm{Pa}]$')
axs1[0][0].set_ylabel(r'Tlak $[\mathrm{Pa}]$')
axs1[0][0].set_title('Čerpání')
#axs1.autoscale(enable=True, axis='x', tight=True)
axs1[1][0].legend()
axs1[0][0].legend()
axs1[1][1].legend()
axs1[0][1].legend()
