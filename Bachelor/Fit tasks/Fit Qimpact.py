def split_fit(X,Y, fitfunction, N=1, maxfev=1000,\
              param_hint=None, split_index=None, fit_smooth=1000, bounds=None, extraplus=0, extraminus=0):
    #split data weather to N pieces, or by split_index and fit each separately 
    #returns array of fit points, and array of fit parameters with standdev
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
                        len(signature(fitfunction).parameters)-1), np.inf)
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
    return(a, np.array([pars, standevs]))