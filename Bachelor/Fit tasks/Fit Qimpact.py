import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

#!!!
def impdata(filename):
    url=f'https://raw.githubusercontent.com/homijan/ML-student-projects/intro-ab/students/bogdaale/gd-profiles/{filename}'
    return url
#!!!

x_Te, Te = np.loadtxt(impdata('Te_gdhohlraum_cm_10ps_TekeV_interp.txt'), usecols=(0, 1), unpack=True)
x_ne, ne = np.loadtxt(impdata('ne_gdhohlraum_cm_ne1e20cm3_interp.txt'), usecols=(0, 1), unpack=True)
x_Zbar, Zbar = np.loadtxt(impdata('Zbar_gdhohlraum_cm_Z_interp.txt'), usecols=(0, 1), unpack=True)

x_Qloc, Qloc = np.loadtxt(impdata('Q_gdhohlraum_microns_10ps_LocalWcm2.txt'), usecols=(0, 1), unpack=True)
x_Qimpact, Qimpact = np.loadtxt(impdata('Q_gdhohlraum_microns_10ps_IMPACTWcm2.txt'), usecols=(0, 1), unpack=True)
x_Qsnb, Qsnb = np.loadtxt(impdata('Q_gdhohlraum_microns_10ps_separatedsnbWcm2.txt'), usecols=(0, 1), unpack=True)

x_Qc7bBGK, Qc7bBGK, Knx = np.loadtxt(impdata('Q_gdhohlraum_cm_10ps_c7b-bgk-Wcm2-clogCHIC.txt'), comments='#', delimiter=',', usecols=(0, 8, 6), unpack=True)
x_Qc7bAWBS, Qc7bAWBS = np.loadtxt(impdata('Q_gdhohlraum_cm_10ps_c7b-awbs-Wcm2-clogCHIC.txt'), comments='#', delimiter=',', usecols=(0, 8), unpack=True)


#changing units um->cm
x_Qloc/=1e4
x_Qimpact/=1e4
x_Qsnb/=1e4


#calculating Te gradient
gradTe=np.gradient(Te, x_Te)


#!!! 3-th task
def fitQloc(X, k):
    #fit function for Qloc profile
    Z, T, gradT = X
    q = -(k/Z)*((Z+0.24)/(Z+4.2))*T**2.5*gradT
    return q
par3, cov3 = curve_fit(fitQloc, (Zbar, Te, gradTe), Qloc,  maxfev = 1000)
standev3=np.sqrt(np.diag(cov3))
print(f'Constant from Qloc profile k = {par3[0]:.1e} ± {standev3[0]:.1e}')


#!!! 4-th task
def fitQimpact(X, alphaC, alphaN):
    #fit function for Qloc profile
    Z, T, gradT = X
    q = -(alphaC/Z)*((Z+0.24)/(Z+4.2))*T**alphaN*gradT
    return q

#!!! 5-th task
print('\n evaluate "nonlocal" fitting constants above subintervals \n')
allfits = np.empty([4,2,x_Te.shape[0]]) #will contain all fits for all N in [1,5]
for N in range(1,5):
    print("     For N = ", N, '\n')
    #split the data we use in fit
    split_Z = np.array_split(Zbar, N)
    split_T = np.array_split(Te, N)
    Split_gradT = np.array_split(gradTe, N)
    split_Qimpact = np.array_split(Qimpact, N)    #split the data we want to fit
    split_x = np.array_split(x_Te, N)             #split the x axes to know at which subinterval we currently are
    subfits = np.array([[],[]])                   #contain fit data for each subinterval
    for sub, _ in enumerate(split_x):             #fitting for each subinterval of x_Te
        pars, covs = curve_fit(fitQimpact, (split_Z[sub], split_T[sub], Split_gradT[sub]), \
                               split_Qimpact[sub],  maxfev = 100000)
        standevs = np.sqrt(np.diag(covs))
        print(f'Pars in subinterval x in <{split_x[sub][0]:.3e} ; {split_x[sub][-1]:.3e}>:','\n',\
              f'alpha_C = {pars[0]:.2e} ± {standevs[0]:.2e}','\n',\
              f'alpha_N = {pars[1]:.2e} ± {standevs[1]:.2e}','\n')
            
        subfits = np.concatenate((subfits, np.array([split_x[sub], \
                                  fitQimpact((split_Z[sub], split_T[sub], Split_gradT[sub]),*pars)])), axis = 1)
    allfits[N-1] = subfits

    
#plot stuff
fontsize = 15.5
plt.rcParams.update({'font.size': fontsize})

fig1, axs1 = plt.subplots(1, 1, figsize=(8, 8))
axs1.plot(x_Te, Te, label="Te")
axs1.set_xlabel('cm')
axs1.set_ylabel('eV')
axs1.set_title('Te')
#axs1.autoscale(enable=True, axis='x', tight=True)
axs1.legend()

fig2, axs2 = plt.subplots(1, 1, figsize=(8, 8))
axs2.plot(x_ne, ne, label="ne")
axs2.set_xlabel('cm')
axs2.set_ylabel('1/cm$-3$')
axs2.set_title('ne')
#axs2.autoscale(enable=True, axis='x', tight=True)
axs2.legend()

fig3, axs3 = plt.subplots(1, 1, figsize=(8, 8))
axs3.plot(x_Zbar, Zbar, label="Zbar")
axs3.set_xlabel('cm')
axs3.set_title('Zbar')
#axs3.autoscale(enable=True, axis='x', tight=True)
axs3.legend()

fig4, axs4 = plt.subplots(1, 1, figsize=(8, 8))
axs4.plot(x_Qloc, Qloc, label="Qloc")
axs4.plot(x_Qimpact, Qimpact, label="Qimpact")
axs4.plot(x_Qsnb, Qsnb, label="Qsnb")
#axs4.plot(1e4 * x_Qc7bAWBS, Qc7bAWBS, label="Qc7b-awbs")
axs4.plot(x_Qc7bAWBS, gaussian_filter1d(Qc7bAWBS,3), label="Qc7b-awbs")
#axs4.plot(1e4 * x_Qc7bBGK, Qc7bBGK, label="Qc7b-bgk")
axs4.plot(x_Qc7bBGK, gaussian_filter1d(Qc7bBGK, 3), label="Qc7b-bgk")
axs4.set_xlabel('cm')
axs4.set_ylabel('W/cm$^2$')
axs4.set_title('Q')
#axs4.autoscale(enable=True, axis='x', tight=True)
axs4.legend(loc='upper left')


fig5, axs5 = plt.subplots(1, 1, figsize=(8, 8))
axs5.plot(x_Qc7bBGK, Knx)
axs5.set_xlabel('cm')
axs5.set_ylabel('[-]')
axs5.set_title(label=r"Knudsen number $Kn_{\mathrm{x}}$")
plt.show()


fig6, axs6 = plt.subplots(2, 2, figsize=(16, 8))
for i in range(4):
    axs6[i//2,i%2].plot(allfits[i,0,:],allfits[i,1,:], label =f'fit for N={i+1}')
    axs6[i//2,i%2].plot(x_Qimpact, Qimpact, label="Qimpact", alpha=0.5)
    axs6[i//2,i%2].set_xlabel('cm')
    axs6[i//2,i%2].set_ylabel('W/cm$^2$')
    axs6[i//2,i%2].legend(loc='upper left')





