import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import lmfit

bez_prop=pd.read_csv("C:/Users/aleks/Projects/VKT/olej_rot/bez_prop.csv", header=0)
bez_prop['p1[Pa]']=bez_prop['p1[Pa]'].str.replace(',','.').astype(float)
s_prop=pd.read_csv('C:/Users/aleks/Projects/VKT/olej_rot/1_s_prop - List 1.csv', header=0)

fig1, axs1 = plt.subplots(1, 1, figsize=(3, 3))
axs1.plot(s_prop['t[s]'], s_prop['p2[Pa]'], label="S proplachovanim")
axs1.plot(bez_prop['t[s]'], bez_prop['p1[Pa]'], label="Bez proplachovani")
axs1.set_xlabel('s')
axs1.set_ylabel('Pa')
axs1.set_title('Čerpání')
#axs1.autoscale(enable=True, axis='x', tight=True)
axs1.legend()
plt.yscale("log")


def fitfunc(x, SV):
    p=np.exp(np.log(96e3)-SV*(x))
    return p
par, cov = curve_fit(fitfunc, bez_prop['t[s]'], bez_prop['p1[Pa]'], p0=1e-2,  maxfev = 100000)
par2, cov2 = curve_fit(fitfunc, s_prop['t[s]'], s_prop['p2[Pa]'], p0=1e-2, maxfev = 100000)
standev=np.sqrt(np.diag(cov))
standev2=np.sqrt(np.diag(cov))
print(par,'  -', standev)
print(par2,'  -', standev2)