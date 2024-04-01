#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:04:31 2023

@author: jep160
"""
#quantifies number of spots, normalized per area, per minute [initiation rate]

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
from glob import iglob

T = 5 #duration of experiment in minutes

#this code takes Tracj duration from all files, removes first frame trajectory and all those that were at the last frame
files = ('/your file directory/')
path = iglob(files, recursive=True)


All= []

for file in path:
    
    
    df = pd.read_csv((file)+'/Track statistics.csv')


    F0 = df.TRACK_START.sort_values().unique()[0] #first frame (always zero) in seconds
    q = df.TRACK_START.sort_values().unique()[1]  #Second frame, also acquisition time
    F1 = df.TRACK_STOP.sort_values().unique()[-1]  #last frame, long of recording in seconds

    AreaC = pd.read_csv((file)+'/Results.csv')
    TA = AreaC.iloc[0,1]
    
    df1 = df.set_index(['TRACK_START'])
    df1 = df1[df1.index>0] #remove trajectories that start at the first frame
    df1 = df1.reset_index()
    df1 = df1.set_index(['TRACK_STOP'])
    df1 = df1[df1.index<F1] #remove trajectories that were cut at the last frame
    df1 = df1.reset_index()
    DD1 = len(df1.TRACK_ID)/TA #number of trajectories per area of cell
    DD  = pd.DataFrame([DD1/T])
    DD.to_csv((file)+'/S.csv', sep=',')   
    print(file)
    
    a = pd.read_csv((file)+'/S.csv')   

path1 = files+'S.csv'
all_rec1 = iglob(path1, recursive=True) 
dataframes1 = (pd.read_csv(f).iloc[:,1] for f in all_rec1)
big_dataframe1 = pd.concat(dataframes1, axis=0, join='outer')
big_dataframe1.to_csv('/Users/Jep160/Desktop/Alfa1.csv')


path = iglob(files, recursive=True)
for file in path:
    
    if os.path.exists((file)+'/S.csv'):
        
        os.remove((file)+'/S.csv')
