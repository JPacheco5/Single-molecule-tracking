#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 11:52:24 2024

@author: jep160
"""

#Model that fit a lognormal distributions and density of molecules by using TrackMate files
import os ## allows change directory##

import matplotlib
import matplotlib.pyplot as plt
import random

import pandas as pd
import numpy as np
import itertools
from PIL import Image
from PIL import ImageDraw
import cv2
import math
from scipy.stats.mstats import gmean
from scipy.optimize import curve_fit
import seaborn; seaborn.set_style('whitegrid')
from glob import iglob


name = "P15"

path = r'/Users/Jep160/Desktop/'+ name+'/*/Spots in tracks statistics.csv'
all_rec = iglob(path, recursive=True) 
dataframes = (pd.read_csv(f).iloc[:,[2,7,15]] for f in all_rec)  #15 for Max Intensity
INT2 = pd.concat(dataframes, axis=1, join='outer')
INT = INT2.loc[:,'MAX_INTENSITY']
INT3 = INT2.loc[:,'TRACK_ID']
INT4 = INT2.loc[:,'POSITION_T']


area = r'/Users/Jep160/Desktop/'+name+'/*/Results.csv'
all_rec1 = iglob(area, recursive=True) 
dataframesA = (pd.read_csv(c).loc[:,['Area']] for c in all_rec1)  
AA = pd.concat(dataframesA, axis=1, join='outer')

ek = INT.columns.values
each = np.arange(len(ek))


GM1 = []
SM1 = []
R2 = []
NM = []
for i in each:
    
    df2 = INT.iloc[:,i].dropna()
    df3 = np.array(df2.tolist())
    df4 = df3.reshape((len(df3),1))

    Xi = np.arange(1,500,1)
    Yi, bins = np.histogram(df4,Xi)
    Y = Yi/Yi.sum()
    X = bins[:-1]


    def logN(x, mu1, sigma1):
        
    
        A = (np.exp(-(np.log(x) - mu1)**2 / (2 * sigma1 **2)) / (x * sigma1 * np.sqrt(2 * np.pi)))
      
        return  A #lognormal function

    params, pcov = curve_fit(logN, X,Y, method="trf", bounds=((0,0),(np.inf,np.inf)), p0=(1,1), maxfev=4000)
    gm = (params[0],params[1]) #gmeans and gSD
    

    gm2 = gm[0]
    sm = gm[1]
    GM1.append(gm2)
    SM1.append(sm)

    print(np.exp(params[0]),np.exp(params[1]))

    
    #R^2 estimation
    residuals = Y - logN(X ,params[0], params[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((Y-np.mean(Y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    R2.append(r_squared)
    
    #GMM

    #model = GeneralMixtureModel.from_samples(LogNormalDistribution, n_components=2, X=df4)
                                       
    #print(np.exp(model.distributions[0].parameters[0]))
    #print(np.exp(model.distributions[1].parameters[0]))



    plt.figure(figsize=(5, 5))  #size of graph
    plt.plot(X, Y, 'r', linewidth=2)
    plt.plot(X, logN(X ,params[0], params[1]),'b', linewidth=2) 
    plt.xlim([-5, 200])
    plt.ylim([0, 0.12])
    #plt.plot(j, model.probability(j), 'g-o')


    plt.show()

for i in each:
    df5 = INT3.iloc[:,i].dropna()
    q = INT4.iloc[:,i].max()
    mol = df5.drop_duplicates().shape[0] #total number of molecules
    areaF = AA.iloc[:,i]
    Nh = (mol/areaF)/q # normalize number of molecules per area and per second
    NM.append(Nh)

    
 
    
GMs = np.reshape(np.exp(GM1), (len(ek),1))
SMs = np.reshape((SM1), (len(ek),1))
Num = np.reshape((NM), (len(ek),1))
R = pd.DataFrame(np.reshape((R2), (len(ek),1)),columns=['R2'])
gsm = pd.DataFrame(SMs, columns=['gSD'])
Gmean = pd.DataFrame(GMs, columns=['Gmean1'])
NumM = pd.DataFrame(Num, columns=['Molecules/µm'])

data = pd.concat([NumM, Gmean, gsm, R], axis=1)
data[data > 1000] = np.nan

data.to_csv('/Users/Jep160/Desktop/DModel1.csv', sep=',') 

print(np.log(data.Gmean1.mean()))
print(data.gSD.mean())

print(data.iloc[:,:2])