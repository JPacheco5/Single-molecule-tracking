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
files = ('/Users/file_location/')
path = iglob(files, recursive=True)




for file in path:
    
    where = 0 #start analysis from what second of experiment  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    start = 0.8 #smallest length of trajectory to analyze, in time (seconds) 
    to = 16 #maximum long of trajectory to analyze (clipped trajectories) in steps
    fitX = 4 #number of points to fit (4-15)
    name = 'Diffusion-1-4'  #experiment <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    name1= 'MSDA'
    
    
    #CALCULATES DIFFUSION COEFFICIENT and ALPHA FROM MSD
    x = pd.DataFrame()  

    #Open the file: in green put the name of the file
    df=pd.read_csv((file)+ '/1Spots in tracks statistics.csv')    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
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

    sf = pd.read_csv((file) + '/1Track statistics.csv')        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    sf1 = sf.sort_values(by=['TRACK_DURATION'])
    sf2 = sf1[['TRACK_ID','TRACK_DURATION','TRACK_DISPLACEMENT', 'TRACK_START']].copy()
    sf2 = sf2[sf2.TRACK_START > where]   #from what second start the analysis
    sf2.to_csv((file)+'/Sort.csv', sep=',')  


    ef = pd.read_csv((file)+ '/Sort.csv')
    ef1 = ef.set_index(['TRACK_DURATION'])
    sets = ef1.loc[start:, 'TRACK_ID']   # Set trajectories to analyze from time 'start' to the end <<<<<<<<<<<<<<<<<<<<<
    M = sets.values.tolist()
    
    DC = []
    alf = []
    R2 = []
    MSD = []
    for j in M:
    
        dfk = df2.loc[j] 
        dfT = dfk.iloc[:to] #clip trajectories at<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    


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
      
            
            MS = pd.DataFrame(msd)
            MSD.append(MS)
                
            
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
            
            
                 
            return DC, alf, R2, MSD

    
        DfC = alldisplacement(dfT, t_step, coords=['POSITION_X', 'POSITION_Y'])
        data = (pd.DataFrame(list(DfC[0:3]))).transpose()
        data.rename(columns={0:'DC', 1:'alpha', 2:'R2', 3:'TRACK_ID'}, inplace=True)
        data1 = DfC[-1]
        MeanSD = pd.concat(data1,axis=1)
        
    MeanSD.columns = M    
    time = pd.DataFrame(t)
    MMSD = pd.concat ([time,MeanSD],axis=1, ignore_index=False)
    MMSD.to_csv((file)+ '/%s.csv' % name1, sep=',', index=False)
    
    data['TRACK_ID'] = M    
    data.to_csv((file)+ '/%s.csv' % name, sep=',', index=False)
    os.remove((file) + '/Sort.csv')

    print((file))
    
import os; os.system('say "It is done"')
    











