#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 12:51:39 2022

@author: JonathanPacheco
"""

#This code quantify the localizations inside of domains by frame in 1 channels. Also counts the number 
#of molecules that fall in domain and normalize by the area out and in domain. It requires the area of the cell

import os ## allows change directory##

import matplotlib
import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import itertools
from PIL import Image
from PIL import ImageDraw
from glob import iglob

start = 0 #minimum length of trajectories to consider for analysis
px = round(0.085, 3) #conversion micrometers to pixels   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
L = 250 #  size of frame (pixels) 


files = ('/Users/Jep160/Desktop/AA/*/*/')
path = iglob(files, recursive=True)

for file in path:
    
    #generate the domain mask
    kf = pd.read_csv((file)+'/mask.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_DURATION'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT']].copy()




    dff2 = pd.read_csv((file)+'/2Spots in tracks statistics.csv')

    q1 = np.round(dff2.loc[1,'POSITION_T']-dff2.loc[0,'POSITION_T'],decimals=3)
    q = float('{0:2f}'.format(q1))


    kList2 = []
    
    for h2 in dff2.TRACK_ID.unique():
    
        ListF2 = np.arange(len(dff2[dff2['TRACK_ID']==h2]))
        kList2.extend(ListF2)
    dff2['nFRAMES'] = kList2
    dff2.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)


    dff2['FRAMES'] = dff2.nFRAMES  #Frames from 1 to N

    df22 = dff2[['TRACK_ID2','POSITION_X','POSITION_Y','POSITION_T','FRAME','FRAMES']].copy() 
    df22['X'] = (np.around(df22['POSITION_X'].div(px), decimals=3)).astype(int)
    df22['Y'] = (np.around(df22['POSITION_Y'].div(px), decimals=3)).astype(int)



    bf2 = sort2.set_index(['TRACK_DURATION'])
    setsB = bf2.loc[start:, 'TRACK_ID2']   # Set trajectories to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = setsB.values.tolist()
    dffb = ((df22.set_index(['TRACK_ID2'])).loc[Wb]).reset_index()

    # Merge for Ch2
    M2 = dffb.merge(domain, on=['X','Y'], how='inner', indicator='pair')
    C2 = M2['POSITION_T'].value_counts()   #localizations inside domian
    C2 = C2.reset_index()
    C2 = C2.sort_values(by=['index'],ascending=True)
    C2.rename(columns={'index':'Time2', 'POSITION_T':'CountsCh2'}, inplace=True)
    C2 = C2.set_index(['Time2'])

    R = pd.concat([C2], axis=1)
    inx = pd.DataFrame(index=np.round(np.arange(0,q*(df22.FRAME.max()+1),q), decimals=2))
    results = R.reindex(inx.index).fillna(0)
    #results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations.csv')


    print(file)
    os.remove((file)+'/domain.csv')

import os; os.system('say "It is done"')