#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 02:58:54 2020

@author: JonathanPacheco
"""
from glob import iglob


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
olderr = np.seterr(all='ignore')



from scipy.optimize import curve_fit
files = ('/file directory/')
path = iglob(files, recursive=True)




for file in path:
    

    start = 0.2 #smallest length of trajectory to analyze, in time (seconds) 
    to = 4 #maximum long of trajectory to analyze (clipped trajectories) in steps
    fitX = 4 #number of points to fit (4-15)
    name = 'Diffusion-PTH-bound'  #experiment <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #CALCULATES DIFFUSION COEFFICIENT and ALPHA FROM MSD

    #Open the file: in green put the name of the file
    df=pd.read_csv((file)+ '/2Spots in tracks statistics.csv')    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #FIX FRAMES
    kList = []
    
    for h in df.TRACK_ID.unique():
        ListF = np.arange(len(df[df['TRACK_ID']==h]))
        kList.extend(ListF)
    df['nFRAMES'] = kList
    df.rename(columns={'FRAME':'old_FRAMES', 'nFRAMES':'FRAME'}, inplace=True)
    
    
    #generete a clean dataframe
    df1 = df[['TRACK_ID','POSITION_X','POSITION_Y','POSITION_T']].copy()

    df2 = df1.set_index(['TRACK_ID'])

    q1 = np.round(df1.iloc[1,3]-df1.iloc[0,3],decimals=3)
    q = float('{0:2f}'.format(q1))

    #selection of tracks by duration

    sf = pd.read_csv((file) + '/2Track statistics.csv')        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    sf1 = sf.sort_values(by=['TRACK_DURATION'])
    sf2 = sf1[['TRACK_ID','TRACK_DURATION','TRACK_DISPLACEMENT']].copy()
    sf2.to_csv((file)+'/Sort.csv', sep=',')  


    ef = pd.read_csv((file)+ '/Sort.csv')
    ef1 = ef.set_index(['TRACK_DURATION'])
    sets = ef1.loc[start:, 'TRACK_ID']   # Set trajectories to analyze from time 'start' to the end <<<<<<<<<<<<<<<<<<<<<
    
    ed = pd.read_csv((file)+ '/ResultsCh2.csv')  ################# to graph interactive molecules   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    setsd = (ed.TRACK_ID2).dropna()               #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<################# to graph interactive molecules
                   ################# to graph interactive molecules
    
    tracksT = pd.DataFrame(sets)                                  #$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    tracksINT = pd.DataFrame(setsd)                                 #$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    tracksINT = tracksINT.rename(columns={'TRACK_ID2':'TRACK_ID'})     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    newM = tracksT.merge(tracksINT, how='left', indicator=True)        #$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    newM = newM[(newM['_merge']=='left_only')].copy()                   #$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    newM = newM.drop(columns='_merge').copy()                          #$$$$$$$$$$$$$$$$$$$ to graph non interacting molecules
    newM1 = newM.TRACK_ID     #trajectories that not interact
    
    #Choose the set of data to analyze
    #M = sets.values.tolist()  #all trajectories
    M = setsd.values.tolist()   #Trajectories witn interaction
    #M = newM1.values.tolist()   #Trajectories that no interact
    DC = []
    alf = []
    R2 = []

    for j in M:
    
        dfk = df2.loc[j] 
        dfT = dfk.iloc[:] #clip trajectories at<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    


    #Parameter input
        N = int(len(dfT))
        max_time = np.float64(N*(q))   
        frames = np.float64(max_time/N)
        t_step = q
        t = np.linspace(q, max_time, N) 

        #function to measure MSD (all displacement)
        def alldisplacement(dfT, t_step, coords=['POSITION_X', 'POSITION_Y']):


            tau = t.copy()
            shifts = (np.round(np.divide(tau,t_step),decimals=3)).astype(int)
            msds_sum = np.zeros(shifts.size)
            delta_inv = np.arange(N)
            delta = delta_inv[N-1::-1]
        
        
            for i, shift in enumerate(np.round(shifts,0)):
                diffs = dfT[coords] - dfT[coords].shift(-shift)
                sqdist = np.square(diffs).sum(axis=1)
                msds_sum[i] = sqdist.sum()
                msd1 = np.divide(msds_sum,delta)
                msd = pd.DataFrame(msd1).dropna()
            
            fitL = fitX  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>number of points to fit 
            Y = msd.loc[:fitL,0]
            X = t[0:len(Y)]
       
            def func(X, D, alpha):
                return 4*D*(X**alpha)
                 
            popt, pcov = curve_fit(func, X, Y, maxfev=4000)

                         
            D = popt[0]
            alpha = popt[1]
            DC.append(D)
            alf.append(alpha)
                
            #R^2 estimation
            residuals = Y - func(X,D,alpha)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((Y-np.mean(Y))**2)
            r_squared = 1 - (ss_res / ss_tot)
            R2.append(r_squared)
        
            return DC, alf, R2

    
        DfC = alldisplacement(dfT, t_step, coords=['POSITION_X', 'POSITION_Y'])
        data = (pd.DataFrame(list(DfC))).transpose()
        data.rename(columns={0:'DC', 1:'alpha', 2:'R2', 3:'TRACK_ID'}, inplace=True)
        
    data['TRACK_ID'] = M    
    data.to_csv((file)+ '/%s.csv' % name, sep=',', index=False)
    os.remove((file) + '/Sort.csv')

    print((file))
    
import os; os.system('say "It is done"')
    











