N=4         #Number of domains
fitoption=6 #Choose one of the options of fitting:
            #1 - constant alphaC, alphaN = 2.5
            #2 - monomic alphaC, alphaN = 2.5
            #3 - constant alphaC, constant alphaN
            #4 - monomic alphaC, constant alphaN 
            #5 - Quadratic alphaC, alphaN = 2.5
            #6 - monomic alphaC, monomic alphaN
xmin = -1   #Set range of interest 
xmax = 1

########
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from scipy.interpolate import CubicSpline
import math
import sys
#%matplotlib ipympl
from numba import jit, cuda
#!!!
def import_from_git(filename):
    url=f'https://raw.githubusercontent.com/homijan/ML-student-projects/intro-ab/students/bogdaale/gd-profiles/{filename}'
    return url
#!!!
#############

x_Te, Te_raw = np.loadtxt(import_from_git('Te_gdhohlraum_cm_10ps_TekeV_interp.txt'), usecols=(0, 1), unpack=True)
x_ne, ne_raw = np.loadtxt(import_from_git('ne_gdhohlraum_cm_ne1e20cm3_interp.txt'), usecols=(0, 1), unpack=True)
x_Zbar, Zbar_raw = np.loadtxt(import_from_git('Zbar_gdhohlraum_cm_Z_interp.txt'), usecols=(0, 1), unpack=True)

x_Qloc, Qloc_raw = np.loadtxt(import_from_git('Q_gdhohlraum_microns_10ps_LocalWcm2.txt'), usecols=(0, 1), unpack=True)
x_Qimpact, Qimpact_raw = np.loadtxt(import_from_git('Q_gdhohlraum_microns_10ps_IMPACTWcm2.txt'), usecols=(0, 1), unpack=True)
x_Qsnb, Qsnb_raw = np.loadtxt(import_from_git('Q_gdhohlraum_microns_10ps_separatedsnbWcm2.txt'), usecols=(0, 1), unpack=True)

x_Qc7bBGK, Qc7bBGK_raw, Knx_raw = np.loadtxt(import_from_git('Q_gdhohlraum_cm_10ps_c7b-bgk-Wcm2-clogCHIC.txt'), comments='#',\
                                         delimiter=',', usecols=(0, 8, 6), unpack=True)
x_Qc7bAWBS, Qc7bAWBS_raw = np.loadtxt(import_from_git('Q_gdhohlraum_cm_10ps_c7b-awbs-Wcm2-clogCHIC.txt'), comments='#',\
                                      delimiter=',', usecols=(0, 8), unpack=True)

# changing units um->cm
x_Qloc/=1e4
x_Qimpact/=1e4
x_Qsnb/=1e4

#################

# Use alphaN
opts=pd.DataFrame(index=pd.Index(['useC','useX','useN','useXN','useX2','useXN1']))
opts['values']=np.full(6, True)
if fitoption in range(7):
    opts.iloc[:]=False
    opts.iloc[fitoption-1]=True
else:
    print(f'Unknown option {fitoption}')
    opts['useC']=True
    
    
xref = x_Te[np.logical_and(x_Te > xmin, x_Te < xmax)]


def getsub(f, x, xref):
    f_cs = CubicSpline(x, f)
    return f_cs(xref)

Te = getsub(Te_raw, x_Te, xref)
ne = getsub(ne_raw, x_ne, xref)
Zbar = getsub(Zbar_raw, x_Zbar, xref)
Qloc = getsub(Qloc_raw, x_Qloc, xref)
Qimpact = getsub(Qimpact_raw, x_Qimpact, xref)
Qsnb = getsub(Qsnb_raw, x_Qsnb, xref)
Qc7bBGK = getsub(Qc7bBGK_raw, x_Qc7bBGK, xref)
Knx = getsub(Knx_raw, x_Qc7bBGK, xref)
Qc7bAWBS = getsub(Qc7bAWBS_raw, x_Qc7bAWBS, xref)

    
#calculating Te gradient
gradTe=np.gradient(Te, xref)  
#!!! 3-th task
def fitQloc(X, k):
    #fit function for Qloc profile
    Z, T, gradT = X
    q = -(k/Z)*((Z+0.24)/(Z+4.2))*T**2.5*gradT
    return q
par3, cov3 = curve_fit(fitQloc, (Zbar, Te, gradTe), Qloc,  maxfev = 1000)
standev3=np.sqrt(np.diag(cov3))
kQSH = par3[0]
print(f'Constant from Qloc profile k = {par3[0]:.1e} ± {standev3[0]:.1e}')

##################

opts=pd.DataFrame(index=pd.Index(['useC','useX','useN','useXN','useX2','useXN1']))
opts['values']=np.full(6, True)
if fitoption in range(7):
    opts.iloc[:]=False
    opts.iloc[fitoption-1]=True
else:
    print(f'Unknown option {fitoption}')
    opts['useC']=True
    
def fitQimpact(X, alpha0, alpha1, alpha2, alpha3, alpha4):
    #fit function for Qloc profile
    x, Z, T, gradT = X
    if (opts.loc['useC'][0]):
        q = -(alpha0 * kQSH / Z) * ((Z + 0.24) / (Z + 4.2)) * T**2.5 * gradT 
    else:
        if (opts.loc['useX'][0]):
            q = ( -((alpha0 + alpha1 * x) * kQSH / Z) * ((Z + 0.24)/(Z+4.2)) *T**2.5 * gradT )
        else:
            if (opts.loc['useN'][0]):
                q = ( -(alpha0 * kQSH / Z) * ((Z + 0.24) / (Z + 4.2)) *T**(2.5 * 1.0 / (1.0 + np.exp(alpha3))) * gradT )
            else: 
                if (opts.loc['useXN'][0]):
                    q = ( -((alpha0 + alpha1 * x) * kQSH / Z) * ((Z + 0.24) / (Z + 4.2))* T**(2.5 * 1.0 / (1.0 + np.exp(alpha3))) * gradT )
                else: 
                    if (opts.loc['useX2'][0]):
                        q = ( -((alpha0 + alpha1 * x + alpha2 * x * x) * kQSH / Z) *((Z + 0.24) / (Z + 4.2)) * T**2.5 * gradT )
                    else:
                        if (opts.loc['useXN1'][0]):
                            q = ( -((alpha0 + alpha1 * x) * kQSH / Z) * ((Z + 0.24) / (Z + 4.2)) *\
                                 T**(2.5 * 1.0 / (1.0 + np.exp(alpha3 + alpha4 * x))) * gradT )
                        else:
                            print("Unknown loss function")
    return q
opts


######################################3

# N = 100
# print("     For N = ", N, '\n')
# # split the data we use in fit
# split_Z = np.array_split(Zbar, N)
# split_T = np.array_split(Te, N)
# split_gradT = np.array_split(gradTe, N)
# split_Qimpact = np.array_split(Qimpact, N)    #split the data we want to fit
# split_x = np.array_split(xref, N)             #split the x axes to know at which subinterval we currently are
# subfits = np.array([[],[]])                   #contain fit data for each subinterval
# for sub, _ in enumerate(split_x):             #fitting for each subinterval of xref
#     pars, covs = curve_fit(fitQimpact, (split_x[sub], split_Z[sub],split_T[sub], split_gradT[sub]), split_Qimpact[sub],  maxfev = 100000)
#     standevs = np.sqrt(np.diag(covs))
#     # print(f'Pars in subinterval x in <{split_x[sub][0]:.3e} ; {split_x[sub][-1]:.3e}>:','\n',\
#     #     f'alpha0 = {pars[0]:.2e} +- {standevs[0]:.2e}','\n',\
#     #     f'alpha1 = {pars[1]:.2e} +- {standevs[1]:.2e}','\n',\
#     #     f'alpha2 = {pars[1]:.2e} +- {standevs[1]:.2e}','\n',\
#     #     f'alpha3 = {pars[1]:.2e} +- {standevs[1]:.2e}','\n',\
#     #     f'alpha4 = {pars[1]:.2e} +- {standevs[1]:.2e}','\n')
#     qsub = fitQimpact((split_x[sub], split_Z[sub], split_T[sub],split_gradT[sub]), *pars)
#     subfits = np.concatenate((subfits, np.array([split_x[sub],qsub])), axis = 1)
    
    
###################################3

def getAlphas(x, Z, T, gradT, Qimpact,  width):

    
    rad=0     #finds out how many points fits in radius range
    suma=0
    while suma<=width:
        suma=xref[rad]-xref[0]
        rad+=1
    
    
    slide=pd.DataFrame(columns=['Start','Center', 'End','alphaC','alphaC_std', 'alphaN', 'alphaN_std']) 
    
    for ind, _ in enumerate(x):
        if ind+2*rad>len(x)-1:
            break
        else:
            pars, covs = curve_fit(fitQimpact, (x[ind:ind+2*rad], Z[ind:ind+2*rad], T[ind:ind+2*rad], gradT[ind:ind+2*rad]),\
                                   Qimpact[ind:ind+2*rad],  maxfev = 100000, bounds=[[-300,-300,-300,-300,-300],[300,300,300,300,300]])
            standevs = np.sqrt(np.diag(covs))

            alphaC=(pars[0:3])
            alphaC_std=(standevs[0:3])
            alphaN_std=(standevs[3:])
            alphaN=(pars[3:])

            slide.loc[len(slide.index)]=[x[ind], x[ind+rad], x[ind+2*rad], alphaC, alphaC_std, alphaN, alphaN_std]
            if ind%100==0:
                print(f"We're at {ind}/{len(x)-2*rad}")  
    return(slide)

slides=getAlphas(xref, Zbar, Te, gradTe, Qimpact,0.1)