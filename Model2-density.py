#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:14:33 2024

@author: jep160
"""

#Model that fit a sum of 2 lognormal distributions and density of molecules by using TrackMate files

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

# Generation of 2 population:

areaF = np.square(27.5) #frame size in microns
path = r'/Users/Jep160/Desktop/X1/*/Spots in tracks statistics.csv'
all_rec = iglob(path, recursive=True) 
dataframes = (pd.read_csv(f).iloc[:,[2,7,15]] for f in all_rec)  #15 for Total Intensity
INT2 = pd.concat(dataframes, axis=1, join='outer')
INT = INT2.loc[:,'MAX_INTENSITY']
INT3 = INT2.loc[:,'TRACK_ID']
INT4 = INT2.loc[:,'POSITION_T']


ek = INT.columns.values
each = np.arange(len(ek))


GM1 = []
PM1 = []
R2 = []
SM1 = []
NM = []

for i in each:
    
    df2 = INT.iloc[:,i].dropna()
    df3 = np.array(df2.tolist())
    df4 = df3.reshape((len(df3),1))

    Xi = np.arange(1,501,1)
    Yi, bins = np.histogram(df4,Xi)
    Y = Yi/Yi.sum()
    X = bins[:-1]


    def logN2(x, mu1, mu2, sigma1, sigma2, P1, P2):
        
        
        A = (P1)*(np.exp(-(np.log(x) - mu1)**2 / (2 * sigma1 **2)) / (x * sigma1 * np.sqrt(2 * np.pi)))
        B = (P2)*(np.exp(-(np.log(x) - mu2)**2 / (2 * sigma2 **2)) / (x * sigma2 * np.sqrt(2 * np.pi)))

      
        return  A+B  #lognormal function
    
    mu1L = 3
    mu1H = 7
    mu1po = 5
    mu1sd = 1
    mu2L = 3
    mu2H = 7
    mu2po = 5
    mu2sd = 1

    
    params, pcov = curve_fit(logN2, X,Y, method="trf", bounds=((mu1L,mu2L,0,0,0,0),(mu1H,mu2H,10,10,1,1)), p0=(mu1po,mu2po,mu1sd,mu2sd,0.1,0.1), maxfev=4000)
    print(params[4]+params[5]) #sum of fractions (must be close to 1)
    gm = (params[0],params[1]) #gmeans
    Ps = (params[4],params[5])
    SD = (params[2],params[3])
    mx = zip(gm,Ps,SD)
    mix = np.array(sorted(mx))  #sort from min to max the gmeans and assign those values to fractions
    gm3 = (mix[:,[0,]]).T
    pm3 = (mix[:,[1,]]).T
    sm3 = (mix[:,[2,]]).T
    GM1.append(gm3)
    PM1.append(pm3)
    SM1.append(sm3)

    print(np.exp(params[0]),np.exp(params[1]))
    print(params[4],params[5])
    
    #R^2 estimation
    residuals = Y - logN2(X ,params[0], params[1],params[2], params[3], params[4], params[5])
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
    plt.plot(X, logN2(X ,params[0], params[1],params[2], params[3], params[4], params[5]),'b', linewidth=2) 
    plt.xlim([-5, 150])
    plt.ylim([0, 0.07])
    #plt.plot(j, model.probability(j), 'g-o')
    plt.plot(X, logN2(X ,params[0], params[1],params[2], params[3], params[4], 0),'o',alpha=0.2) 
    plt.plot(X, logN2(X ,params[0], params[1],params[2], params[3], 0, params[5]),'o',alpha=0.2) 


    plt.show()
    
    
for i in each:
    df5 = INT3.iloc[:,i].dropna()
    q = INT4.iloc[:,i].max()
    mol = df5.drop_duplicates().shape[0] #total number of molecules
    Nh = (mol/areaF)/q # normalize number of molecules per area and per second
    NM.append(Nh)
    
GMs = np.reshape(np.exp(GM1), (len(ek),2))
PMs = np.reshape((PM1), (len(ek),2))
SMs = np.reshape((SM1), (len(ek),2))
Num = np.reshape((NM), (len(ek),1))
#R = pd.DataFrame(np.reshape((R2), (len(ek),1)),columns=['R2'])
Frac = pd.DataFrame(PMs, columns=['F1', 'F2'])
Gmean = (pd.DataFrame(GMs, columns=['Gmean1', 'Gmean2'])).round(decimals=0)
gSM = pd.DataFrame(SMs, columns=['gSD1', 'gSD2'])
NumM = pd.DataFrame(Num, columns=['Molecules/µm'])


data = pd.concat([Frac, Gmean, gSM, NumM], axis=1)
data['P1'] = np.where((Gmean['Gmean1'] != Gmean['Gmean2']), data['F1'], 1)
data['P2'] = 1 - data['P1']
data['GM1'] = np.where((Gmean['Gmean1'] != Gmean['Gmean2']), data['Gmean1'], data['Gmean1'])
data['GM2'] = np.where((data['P1'] != 1), data['Gmean2'], np.nan)
data[data > 600] = np.nan



data.to_csv('/Users/Jep160/Desktop/Model2.csv', sep=',')
