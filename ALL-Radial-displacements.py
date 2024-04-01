#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 12:06:31 2021

@author: JonathanPacheco
"""

#This code calculates Radial displacements
#You must update time adquisition in seccion of parameter input and t

#libraries
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
from scipy.optimize import curve_fit


files = ('/Users/your file directory/')
path = iglob(files, recursive=True)

for file in path:
    
    #parameters
    where = 0 #start analysis from what second of experiment
    until = 300 #last second of experiment to analyze
    start = 0.2  #smallest length of trajectory to analyze, in time (seconds) 
    #to = 16 #maximum length of trajectory to analyze (clipped trajectories) in steps        

    #Open the file: in green put the name of the file
    df=pd.read_csv((file)+'/1Spots in tracks statistics.csv')


    #generete a clean dataframe
    df1 = df[['TRACK_ID','POSITION_X','POSITION_Y','POSITION_T']].copy()
    df2 = df1.set_index(['TRACK_ID'])
    q1 = np.round(df1.iloc[1,3]-df1.iloc[0,3],decimals=3)
    q = float('{0:2f}'.format(q1))


    #selection of tracks by duration

    sf = pd.read_csv((file) + '/1Track statistics.csv')
    sf1 = sf.sort_values(by=['TRACK_DURATION'])
    sf2 = sf1[['TRACK_ID','TRACK_DURATION','TRACK_DISPLACEMENT', 'TRACK_START']].copy()
    sf2 = sf2[sf2.TRACK_START > where]
    sf3 = sf2[sf2.TRACK_START < until]
    sf3.to_csv((file)+'/Sort.csv', sep=',')  
   

    ef = pd.read_csv((file)+ '/Sort.csv')
    ef1 = ef.set_index(['TRACK_DURATION'])
    sets = ef1.loc[start:, 'TRACK_ID']   # Set trajectories to analyze from time X to the end <<<<<<<<<<<<<<<<<<<<<
    M = sets.values.tolist()
    RD = []



    for j in M:
    
        dfk = df2.loc[j] 
        dfT = dfk.iloc[:] #  enter number of points to analyze.  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        N = int(len(dfT))
        max_time = np.float64(N*(q))  
        frames = np.float64(max_time/N)
        t_step= frames

        data = pd.DataFrame({'N':[N],'max_time':[max_time],'frames':[frames]})

        t= np.round(np.linspace(q, max_time, N), decimals=3)



    #function to get radial displacement (all displacement)
        def radial(dfT, coords=['POSITION_X', 'POSITION_Y']):


            tau = t.copy()
            shifts = (np.round(np.divide(tau,t_step),decimals=3)).astype(int)
            radials = []

            for i, shift in enumerate(np.round(shifts,0)):
                diffs = np.array(dfT[coords] - dfT[coords].shift(-shift))
                sqdist = np.square(diffs).sum(axis=1)
                r = np.sqrt(sqdist)
                radials.append(r)


    
            radial_disp = pd.DataFrame({'radials':radials})
            return radials
    
    
        radial_d = radial(dfT, coords=['POSITION_X', 'POSITION_Y'])

        radd = pd.DataFrame.from_records(radial_d) #horizontal
        rad = radd.transpose() #vertical

        b = rad
        RD.append(b)
        b.to_csv((file)+ '/radial_dis.csv', sep=',',mode='a', index=True)
           
    JP1 = pd.concat(RD)
    JP2 = JP1.iloc[:,0:4]
    JP2.columns = t[0:4]

    JP2.to_csv((file)+ '/radial resultsR.csv', sep=',',mode='a', index=True)
    os.remove((file)+ '/Sort.csv')
    os.remove((file)+ '/radial_dis.csv')
    print((file))
   

import os; os.system('say "It is done"')
print('++++++++++++++++++++++++++++++++++++++++++Radial displacement analysis has finished++++++++++++++++++++++++++++++++++++++++++')