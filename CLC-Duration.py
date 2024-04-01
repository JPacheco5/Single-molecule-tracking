#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:23:56 2022

@author: JonathanPacheco
"""
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

#this code takes Track duration from all files, removes first frame trajectory and all those that were at the last frame
files = ('/Your file directory/')
path = iglob(files, recursive=True)


All= []

for file in path:
    
    
    df = pd.read_csv((file)+'/Track statistics.csv')


    F0 = df.TRACK_START.sort_values().unique()[0] #first frame (always zero) in seconds
    q = df.TRACK_START.sort_values().unique()[1]  #Second frame, also acquisition time
    F1 = df.TRACK_STOP.sort_values().unique()[-1]  #last frame, long of recording in seconds


    df1 = df.set_index(['TRACK_START'])
    df1 = df1[df1.index>0] #remove trajectories that start at the first frame
    df1 = df1.reset_index()
    df1 = df1.set_index(['TRACK_STOP'])
    df1 = df1[df1.index<F1] #remove trajectories that were cut at the last frame
    df1 = df1.reset_index()
    DD = df1.TRACK_DURATION
    All.append(DD)
    print(file)
    
DT = pd.concat(All, axis=1)
DT.to_csv('/your file directory/CLC_duration.csv', sep=',')    

    


    