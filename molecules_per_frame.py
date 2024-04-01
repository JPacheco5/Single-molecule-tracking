#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 13:12:36 2022

@author: jep160
"""

#This code quantify the number of molecules or trajectories normalized by seconds and by area of cell

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

files = ('/Your file directory/')
path = iglob(files, recursive=True)

for file in path:
    
    AreaC = pd.read_csv((file)+'/Results.csv')  
    A = AreaC.iloc[0,1]               # Area of cell (microns)
    
    
    
    #files:
    pf = pd.read_csv((file)+'/2Spots in tracks statistics.csv')
    bf = pd.read_csv((file)+'/2Track statistics.csv')
    
    time = pf.POSITION_T.max()  #seconds that last the experiment
    
    #calculation by trajectories
    tracks = bf.TRACK_ID.max() #number of trajectories
    T_A = tracks/A   #normalization by area
    T_A_norm = T_A/time  #normalization per second
    
    
    #calculation by detections
    mol = pf.TRACK_ID.count()  # total number of molecules
    M_F = mol/pf.POSITION_T.max() # Mean number of molecules per second
    M_F_norm = M_F/A           # Mean number of molecules per frame normalized by area of cell
    print(M_F_norm)


