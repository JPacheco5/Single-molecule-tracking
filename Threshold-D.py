#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:39:29 2024

@author: jep160
"""
#correction diffussion values from MSD using threshold TH

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


#This code remove displacements based on control of immobility

files = ('/Your directory/')
path = iglob(files, recursive=True)



for file in path:
    
    #parameters
    TH = [0.00082]
    

    #Open the file: in green put the name of the file
    df=pd.read_csv((file)+'/Diffusion-1-4.csv')
    
    
    df1 = df.iloc[:,:][df.iloc[:, 0] > TH[0]]

    
    new = pd.concat([df1], axis=1)

    new = new.reset_index(drop=True)
    new.to_csv((file)+'/Diffusion-1-4p.csv')

