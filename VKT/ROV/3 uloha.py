# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 23:26:15 2022

@author: aleks
"""

#zpracovani ulohy z praktik. Mam spektrum zareni cesia-137, hledam fotopik a fituju gausem


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib.ticker import AutoMinorLocator

#format dat: skut hodnota;namerena hodnota;chyba skut;chyba namer
path='C:/Users/aleks/Projects/arameic-mishmash/VKT/olej_rot/'
dataset=pd.read_csv(f'{path}Sef.csv', header=0, delimiter=';', decimal=',')
dataset=dataset.sort_values('p [Pa]')



def fce(x,a,b):
    y = (98e3/x[0])*4.75e-2*a/x[1]+b
    return y
c=20
yaxes=dataset['S_ef'][:];xaxes=dataset['p [Pa]'][:]
par, cov = curve_fit(fce, dataset[['p [Pa]','tau']].values.T, yaxes, maxfev = 1000,\
                     p0=[-15,0], bounds=[[-20,-np.inf],[-0.3,np.inf]]) #!!!
xaxes=np.linspace(dataset.iloc[0]['p [Pa]']-c/6,dataset.iloc[-1]['p [Pa]']+c,100) #aby fitovaci krivka vypadala lip

aaxes=np.linspace(dataset.iloc[0]['tau']-c/6,dataset.iloc[-1]['tau']+c,100)
#pars, covs = curve_fit(fitQimpact, (split_Z[sub], split_T[sub], Split_gradT[sub]), \
#                       split_Qimpact[sub],  maxfev = 100000)



standev=np.sqrt(np.diag(cov))

def plot_param():       #чудо чудное
    FONT = 'Times New Roman'   
    plt.ylabel(r'Sef', fontsize = 15)
    plt.xlabel(r'p[Pa]',  fontsize = 15)
    # plt.xlim(0.,0.5)   #!!!!!
    # plt.ylim(0.,0.5)
    #graphic parameters
    plt.xticks(fontname=FONT, fontweight = 'bold', fontsize = 12)
    plt.yticks(fontname=FONT, fontweight = 'bold', fontsize = 12)
    
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

#vytvarim automaticke popisky
label='$'+str(round(par[0],2))+'x'+ str(round(par[1],2))+'$'

fig, ax=plt.subplots()

fit=ax.plot(xaxes,fce((xaxes,aaxes),*par),color='red', linewidth=2,\
            label=r'$y=(98800/p)(4.75a/\tau)+b$')


ax.legend()

plt.tight_layout()


dataline=ax.scatter(dataset['p [Pa]'], dataset['S_ef'][:],marker='+',\
                    color='black', alpha=0.6, s=40, )

print("Pro první pik")
for i in range(len(par)):
    print(i+1,'. Parameter =', round(par[i], 4), '\t with err ', round(standev[i], 4))

plot_param()   