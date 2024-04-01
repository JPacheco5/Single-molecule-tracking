#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 17:41:02 2022

@author: jep160
"""

#This code quantify the localizations inside of domains by frame in 1 channels. It uses Stacks with increments for longer time-lapses. Also counts the number 
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
#L = 250 #  size of frame (pixels) 

files = ('/your directory/')
path = iglob(files, recursive=True)

for file in path:
    
    #generate the domain mask1
    kf = pd.read_csv((file)+'/mask1.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 0) & (sort2['TRACK_STOP'] <= 30 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations1.csv')
    os.remove((file)+'/domain1.csv')


    #generate the domain mask2  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask2.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 30) & (sort2['TRACK_STOP'] <= 60 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations2.csv')
    os.remove((file)+'/domain1.csv')


    #generate the domain mask3  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask3.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 60) & (sort2['TRACK_STOP'] <= 90 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations3.csv')
    os.remove((file)+'/domain1.csv')


    #generate the domain mask4  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask4.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 90) & (sort2['TRACK_STOP'] <= 120 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations4.csv')
    os.remove((file)+'/domain1.csv')


    #generate the domain mask5  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask5.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 120) & (sort2['TRACK_STOP'] <= 150 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations5.csv')
    os.remove((file)+'/domain1.csv')


    #generate the domain mask6  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask6.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 150) & (sort2['TRACK_STOP'] <= 180 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations6.csv')
    os.remove((file)+'/domain1.csv')
    
    
    #generate the domain mask7  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask7.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 180) & (sort2['TRACK_STOP'] <= 210 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations7.csv')
    os.remove((file)+'/domain1.csv')
    
    #generate the domain mask8  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask8.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 210) & (sort2['TRACK_STOP'] <= 240 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations8.csv')
    os.remove((file)+'/domain1.csv')



    #generate the domain mask9  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask9.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 240) & (sort2['TRACK_STOP'] <= 270 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT

    inL.to_csv((file)+'/Localizations9.csv')
    os.remove((file)+'/domain1.csv')
    
    
    #generate the domain mask10  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    kf = pd.read_csv((file)+'/mask10.txt', sep='\t', engine='python', header=None, skiprows=1)
    kf1 = kf.astype(bool)
    kf2 = kf1.values
    kf3 = np.argwhere(kf2)
    kf4 = pd.DataFrame(kf3)
    kf5 = kf4.to_csv((file)+'/domain1.csv',header=['Y','X'])
    domain = pd.read_csv((file)+'/domain1.csv')
    AreaC = pd.read_csv((file)+'/Results.csv')  
    Total_A = AreaC.iloc[0,1]                    # Total number of pixels in domain image converted to area in micrometers
    Domain_A = (len(domain['X']))*(px**2)               # Number of domain pixels converted to area in micrometers
    nDomain_A = Total_A - Domain_A                  # Number of non domain pixels converted to area in micrometers



    bf = pd.read_csv((file)+'/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_STOP'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT','TRACK_STOP','TRACK_START']].copy()




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



    bf2 = sort2.loc[(sort2['TRACK_START'] >= 270) & (sort2['TRACK_STOP'] <= 296 )]  # Set start-Stop time frame to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = bf2.TRACK_ID2.unique()
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
    results.to_csv((file)+'/Loc-frame.csv')
    
    D2i = M2.TRACK_ID2.value_counts().sum()#total number of molecules in domain Ch2    
     
    D2 = D2i/Domain_A
    
    D2o = len(dffb)-D2i
    D22 = D2o/nDomain_A #total number of molecules out domain Ch2
    
    
    
    inL = pd.DataFrame({'molIN':[D2i],'Ch2IN':[D2],'molOUT':[D2o],'Ch2OUT':[D22]}) 
    inL['ratio1-abs'] = inL.molIN/inL.molOUT  # ratio of detection inside/outside
    inL['ratio2-norm'] = inL.Ch2IN/inL.Ch2OUT  #ratio of normalized detections (area) inside/outside

    inL.to_csv((file)+'/Localizations10.csv')
    os.remove((file)+'/domain1.csv')
        
        
        

    M1 = pd.read_csv((file)+'/Localizations1.csv')
    M2 = pd.read_csv((file)+'/Localizations2.csv')
    M3 = pd.read_csv((file)+'/Localizations3.csv')
    M4 = pd.read_csv((file)+'/Localizations4.csv')
    M5 = pd.read_csv((file)+'/Localizations5.csv')
    M6 = pd.read_csv((file)+'/Localizations6.csv')
    M7 = pd.read_csv((file)+'/Localizations7.csv')
    M8 = pd.read_csv((file)+'/Localizations8.csv')
    M9 = pd.read_csv((file)+'/Localizations9.csv')
    M10 = pd.read_csv((file)+'/Localizations10.csv')
    
    Fin = pd.concat ([M1,M2,M3,M4,M5,M6,M7,M8,M9,M10])
    Fin.to_csv((file)+'/Localizations_intervals.csv')
    
    os.remove((file)+'/Localizations1.csv')
    os.remove((file)+'/Localizations2.csv')
    os.remove((file)+'/Localizations3.csv')
    os.remove((file)+'/Localizations4.csv')
    os.remove((file)+'/Localizations5.csv')
    os.remove((file)+'/Localizations6.csv')
    os.remove((file)+'/Localizations7.csv')
    os.remove((file)+'/Localizations8.csv')
    os.remove((file)+'/Localizations9.csv')
    os.remove((file)+'/Localizations10.csv')
    
    
    print(file)
    