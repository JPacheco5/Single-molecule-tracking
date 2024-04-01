#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 20:27:57 2022

@author: JonathanPacheco
"""
#This code quantifies the time that 2 molecules are colocalizing in a determinated distance (Res)/ also measure the time that 
#each molecule moves freely before interaction and the ratio of molecules that scape from interaction
#This code measure the total time of interaction and the duration of transient interactions, as well the number of transient interactions
#check precision of measurement (line 24)

import os ## allows change directory##
from glob import iglob
import matplotlib
import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import itertools
from PIL import Image
from PIL import ImageDraw

Res = 0.08  #resolution of colocalization in microns (diffusion at one time frame or FWHM)
start = 1.5 #minimum length of trajectories to consider for analysis
px = round(0.085, 3) #conversion micrometers to pixels   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
scale = 100 #usually we use a value of 100
L = 100 * scale # multiply the size of frame (pixels) by the scale
pk = 3000 #final size of image

files = ('/your directory/')
path = iglob(files, recursive=True)

for file in path:
    

    af = pd.read_csv((file)+ '/1Track statistics.csv')
    af1 = af.sort_values(by=['TRACK_DURATION'])
    af1.rename(columns={'TRACK_ID':'TRACK_ID1'}, inplace=True)
    sort1 = af1[['TRACK_ID1','TRACK_DURATION','TRACK_DISPLACEMENT']].copy()


    bf = pd.read_csv((file)+ '/2Track statistics.csv')
    bf1 = bf.sort_values(by=['TRACK_DURATION'])
    bf1.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)
    sort2 = bf1[['TRACK_ID2','TRACK_DURATION','TRACK_DISPLACEMENT']].copy()



    dff1 = pd.read_csv((file)+ '/1Spots in tracks statistics.csv')
    dff2 = pd.read_csv((file)+ '/2Spots in tracks statistics.csv')

    q1 = np.round(dff1.loc[1,'POSITION_T']-dff1.loc[0,'POSITION_T'],decimals=3)
    q = float('{0:2f}'.format(q1))


    kList1 = []
    
    for h1 in dff1.TRACK_ID.unique():
    
        ListF1 = np.arange(len(dff1[dff1['TRACK_ID']==h1]))
        kList1.extend(ListF1)
    dff1['nFRAMES'] = kList1
    dff1.rename(columns={'TRACK_ID':'TRACK_ID1'}, inplace=True)

    kList2 = []
    
    for h2 in dff2.TRACK_ID.unique():
    
        ListF2 = np.arange(len(dff2[dff2['TRACK_ID']==h2]))
        kList2.extend(ListF2)
    dff2['nFRAMES'] = kList2
    dff2.rename(columns={'TRACK_ID':'TRACK_ID2'}, inplace=True)



    dff1['FRAMES'] = dff1.nFRAMES  #Frames from 1 to N
    dff2['FRAMES'] = dff2.nFRAMES  #Frames from 1 to N

    df11 = dff1[['TRACK_ID1','POSITION_X','POSITION_Y','POSITION_T','FRAME','FRAMES']].copy() 
    df22 = dff2[['TRACK_ID2','POSITION_X','POSITION_Y','POSITION_T','FRAME','FRAMES']].copy() 


    af2 = sort1.set_index(['TRACK_DURATION'])
    setsA = af2.loc[start:, 'TRACK_ID1']   # Set trajectories to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wa = setsA.values.tolist()
    dffa = ((df11.set_index(['TRACK_ID1'])).loc[Wa]).reset_index()


    bf2 = sort2.set_index(['TRACK_DURATION'])
    setsB = bf2.loc[start:, 'TRACK_ID2']   # Set trajectories to analyze<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    Wb = setsB.values.tolist()
    dffb = ((df22.set_index(['TRACK_ID2'])).loc[Wb]).reset_index()



    df1 = dffa.sort_values(['FRAME']).set_index('FRAME')
    df2 = dffb.sort_values(['FRAME']).set_index('FRAME')


    fr1 = pd.DataFrame(list(dffa.FRAME.unique())) #it uses a list of Frames from df1
    fr2 = pd.DataFrame(list(dffb.FRAME.unique()))
    fr1.columns=['fr']
    fr2.columns=['fr']


    fr = fr1.merge(fr2, on=['fr'], how='inner', indicator='pair') #merge between frames of each channel
    coloc = []

    for f in fr.fr:
    
        t1 = df1.loc[f]
        t2 = df2.loc[f]

        if t1.size == 5 or t2.size == 5:  #this line avoid empty tables
            continue


        dfs = (t1.assign(T=1).merge(t2.assign(T=1), on='T', suffixes=("","_t2")) #prepare combination of all positions to calculate (T1 column)
        .assign(xdiff=lambda x: x['POSITION_X']-x['POSITION_X_t2'])
        .assign(ydiff=lambda x: x['POSITION_Y']-x['POSITION_Y_t2'])
        .query("xdiff>=-@Res and xdiff<=+@Res and ydiff>=-@Res and ydiff<=+@Res")  #Precision of measurement
        .drop(columns=["T"])
        )
        coloc.append(dfs)
    
    
    P_coloc = pd.concat(coloc)
    Dat = P_coloc.sort_values(['TRACK_ID1','FRAMES'], ascending=[True, True]) #sort Tracks_ID1 and frames
    Datp = P_coloc.sort_values(['TRACK_ID2','FRAMES_t2'], ascending=[True, True]) #sort Tracks_ID2 and frames

    count1 = pd.DataFrame()
    count2 = pd.DataFrame()

    Dat1 = Dat.set_index(['TRACK_ID1'])
    Dat2 = Datp.set_index(['TRACK_ID2'])


    #interacting time 
    for D1 in Dat.TRACK_ID1.unique():
        if Dat1.loc[D1].FRAMES.size == 1:
            continue
        
        CQ = []
        C1p = 1
        ms1 = Dat1.loc[D1].reset_index().FRAMES
        for i in range(len(ms1)):
            if i != (len(ms1) - 1):
                diff = ms1[i+1] - ms1[i]
                if diff <=3:   #here it is established a threshold of lossing by X frames the localization
                    C1p += diff
                else:
                    CQ.append(C1p)
                    C1p = 1
            else:
                CQ.append(C1p)
        
        Dp1 = pd.DataFrame(CQ)*q
        Dp1['TRACK_ID1'] = D1
        count1 = count1.append(Dp1)


    
    for D2 in Datp.TRACK_ID2.unique():
        if Dat2.loc[D2].FRAMES_t2.size == 1:
            continue

        CQ2 = []
        C2p = 1
        ms2 = Dat2.loc[D2].reset_index().FRAMES_t2
        for i2 in range(len(ms2)):
            if i2 != (len(ms2) - 1):
                diff2 = ms2[i2+1] - ms2[i2]
                if diff2 <=3:    #here it is established a threshold of lossing by 1 frame the localization
                    C2p += diff2
                else:
                    CQ2.append(C2p)
                    C2p = 1
            else:
                CQ2.append(C2p)
        
        Dp2 = pd.DataFrame(CQ2)*q
        Dp2['TRACK_ID2'] = D2
        count2 = count2.append(Dp2)
        
        

    count1.rename(columns={0:'Int_time'}, inplace=True)
    count2.rename(columns={0:'Int_time'}, inplace=True)
    if count1.size == 0:
        continue

    if count2.size== 0:
        continue

    count1 = count1.reset_index(drop=True)
    count2 = count2.reset_index(drop=True)

    Dat11 = Dat1.loc[count1.TRACK_ID1.unique()].reset_index()
    Dat22 = Dat2.loc[count2.TRACK_ID2.unique()].reset_index()

    #time pre-interaction
    CH1 = q*(Dat11.groupby('TRACK_ID1').first().reset_index().FRAMES ) #frames at molecules Ch1 interact with CH2
    CH2 = q*(Dat22.groupby('TRACK_ID2').first().reset_index().FRAMES_t2) #frames at molecules Ch1 interact with CH2
    C1 = Dat11.groupby('TRACK_ID1').first().reset_index().TRACK_ID1 
    C2 = Dat22.groupby('TRACK_ID2').first().reset_index().TRACK_ID2 
    ch1 = pd.concat([CH1,C1], axis=1)
    ch2 = pd.concat([CH2,C2], axis=1)
    ch1 = ch1.reset_index()
    ch2 = ch2.reset_index()

    #ratio of molecules that interact from the total
    ratio_ch1 = [len(Dat.TRACK_ID1.unique()) / len(dffa.TRACK_ID1.unique())]  #ratio of molecules that interact from total molecules
    ratio1 = pd.DataFrame(ratio_ch1, columns=['Fraction_int'])

    ratio_ch2 = [len(Dat.TRACK_ID2.unique()) / len(dffb.TRACK_ID2.unique())]  #ratio of molecules that interact from total molecules
    ratio2 = pd.DataFrame(ratio_ch2, columns=['Fraction_int'])

    #Scape molecules: calculate the ratio of molecules that scape from interaction (Channel 2)

    Ch2m = Datp[['TRACK_ID2','POSITION_X_t2','POSITION_Y_t2','POSITION_T_t2', 'FRAMES_t2']].copy() #dataframe containing interacting segments of trajectories
    Ch2m.rename(columns={'POSITION_X_t2':'POSITION_X','POSITION_Y_t2':'POSITION_Y', 'POSITION_T_t2':'POSITION_T'  }, inplace=True)
    dffk = ((dffb.set_index(['TRACK_ID2'])).loc[count2.TRACK_ID2.unique()]).reset_index()  #dataframe containing interacting trajectories (whole)
    nD = dffk.merge(Ch2m, on=['TRACK_ID2','POSITION_X', 'POSITION_Y', 'POSITION_T'], how='left', indicator='pair')
    nD1 = nD.set_index('TRACK_ID2')
    SCP = []

    for EA in nD.TRACK_ID2.unique():
    
    
        DD = nD1.loc[EA].reset_index()
        sel = DD.notna()[::-1].idxmax()
        scape = q*(sel.loc['FRAMES']-sel.loc['FRAMES_t2'])
        if scape > 3*q:  # three times the frame exposure
            scape = 1
            SCP.append(scape)
        else:
            scape = 0
            SCP.append(scape)
        

    release = pd.DataFrame(SCP)
    ratioS = np.divide(release.sum(),release.count())
    Diss = pd.DataFrame(ratioS, columns=['Diss-ratioCH2'])




    #Scape molecules: calculate the ratio of molecules that scape from interaction (Channel 1)

    Ch1m = Dat[['TRACK_ID1','POSITION_X','POSITION_Y','POSITION_T', 'FRAMES']].copy() #dataframe containing interacting segments of trajectories
    
    dffk1 = ((dffa.set_index(['TRACK_ID1'])).loc[count1.TRACK_ID1.unique()]).reset_index()  #dataframe containing interacting trajectories (whole)
    nDc = dffk1.merge(Ch1m, on=['TRACK_ID1','POSITION_X', 'POSITION_Y', 'POSITION_T'], how='left', indicator='pair')
    nDc1 = nDc.set_index('TRACK_ID1')
    SCP1 = []

    for EA1 in nDc.TRACK_ID1.unique():


        DD1= nDc1.loc[EA1].reset_index()
        sel1 = DD1.notna()[::-1].idxmax()
        scape1 = q*(sel1.loc['FRAMES_x']-sel1.loc['FRAMES_y'])
        if scape1 > 3*q:  # three times the frame exposure
            scape1 = 1
            SCP1.append(scape1)
        else:
            scape1 = 0
            SCP1.append(scape1)
    

    release1 = pd.DataFrame(SCP1)
    ratioS1 = np.divide(release1.sum(),release1.count())
    Diss1 = pd.DataFrame(ratioS1, columns=['Diss-ratioCH1'])


    #number of interactions by molecule
    count1['Num_int'] = count1.groupby('TRACK_ID1').count().reset_index().Int_time
    count2['Num_int'] = count2.groupby('TRACK_ID2').count().reset_index().Int_time

    #sum of total time of interaction per trajectory
    TT1 = (count1.groupby('TRACK_ID1').sum().Int_time).reset_index(drop=True)
    TT2 = (count2.groupby('TRACK_ID2').sum().Int_time).reset_index(drop=True)
    TT1 = pd.DataFrame(TT1)
    TT2 = pd.DataFrame(TT2)
    TT1 = TT1.rename(columns={'Int_time':'Total_time'})
    TT2 = TT2.rename(columns={'Int_time':'Total_time'})
    
    #total number of molecules
    Mol1 = pd.DataFrame([len(dffa.TRACK_ID1)])
    Mol1 = Mol1.rename(columns={0:'#Mol1'})
    Mol2 = pd.DataFrame([len(dffb.TRACK_ID2)])
    Mol2 = Mol2.rename(columns={0:'#Mol2'})
    
    
    Result1 = pd.concat([ch1,count1, TT1, ratio1, Diss1, Mol1], axis=1)  
    Result1 = Result1.drop('index',1)
    Result2 = pd.concat([ch2,count2,TT2, ratio2,Diss, Mol2], axis=1)
    Result2 = Result2.drop('index',1)
    Result1.rename(columns={'FRAMES':'pre_int'}, inplace=True)
    Result2.rename(columns={'FRAMES_t2':'pre_int'}, inplace=True)


    Result1.to_csv((file)+'/ResultsCh1.csv', sep=',')  #time molecules CH1 before interaction with CH2
    Result2.to_csv((file)+'/ResultsCh2.csv', sep=',')  #time molecules CH2 before interaction with CH1
    Dat.to_csv((file)+'/P_coloc.csv', sep=',') #Dataframe with sorted tracks_ids colocalizing    
    print(file)

#----------------------print

    Tracks1 = Dat.TRACK_ID1.unique()  #FULL TRAJECTORIES
    Tracks2 = Dat.TRACK_ID2.unique()  #FULL TRAJECTORIES



    #Dat['seg1'] = (Dat['TRACK_ID1']+(Dat['ct']/1000)).astype(float) #Generates new TRACKS ID by summing the frames divided by 10000 to separete non consecutive frames
    #Dat1 = Dat[Dat.seg1.isin(Dat.seg1.value_counts()[Dat.seg1.value_counts()>=2].index)] #is taking at least X steps

    #Dat['ct2'] = Dat.FRAMES_t2.diff().ne(1).cumsum()
    #Dat['seg2'] = (Dat['TRACK_ID2']+(Dat['ct2']/1000)).astype(float) #Generates new TRACKS ID by summing the frames divided by 10000 to separete non consecutive frames
    #Dat2 = Dat[Dat.seg2.isin(Dat.seg2.value_counts()[Dat.seg2.value_counts()>=2].index)] #is taking at least X steps



    #Tracks1 = Dat1.seg1.unique()    #PARTIAL TRAJECTORIES
    #Tracks2 = Dat2.seg2.unique()    #PARTIAL TRAJECTORIES


    print1 = (df11.set_index(['TRACK_ID1'])).loc[Tracks1] #FULL TRAJECTORIES
    print2 = (df22.set_index(['TRACK_ID2'])).loc[Tracks2] #FULL TRAJECTORIES


    #print1 = (Dat1.set_index(['seg1'])).loc[Tracks1] #PARTIAL TRAJECTORIES
    #print2 = (Dat2.set_index(['seg2'])).loc[Tracks2] #PARTIAL TRAJECTORIES



    print1['Xi'] = np.round(print1['POSITION_X'].div(px), decimals=3)*(scale)
    print1['Yi'] = np.round(print1['POSITION_Y'].div(px), decimals=3)*(scale)

    tracks_1 = print1.groupby('TRACK_ID1')

    img = Image.new("RGB", (L,L), "black")   # IMPORTANT to set dimensions of new image <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    draw = ImageDraw.Draw(img)
    #colors = itertools.cycle(('green','green'))
    #"Red", "Chartreuse", "DarkOrange", "DarkMagenta", "gold",'gray','navy','plum','deepskyblue','purple','maroon'
    for e, group in tracks_1:
        coords = tuple(zip(group.Xi, group.Yi))
        #draw.line(coords, width=25, fill=next(colors))  #to draw each trajectory with different colours
        draw.line(coords,width=25, fill='green')    #to draw each trajectory with the same colour
        draw.ellipse([(coords[0]),(tuple([50+x for x in coords[0]]))], fill='deepskyblue')


    maxsize = (pk, pk)
    img = img.resize(maxsize, Image.BICUBIC)
    #img.show()
    traj = img.save((file)+ "/TRACKS1.tif","tiff", transparency=0)


    print2['Xi'] = np.round(print2['POSITION_X'].div(px), decimals=3)*(scale)   #FULL TRAJECTORIES
    print2['Yi'] = np.round(print2['POSITION_Y'].div(px), decimals=3)*(scale)   #FULL TRAJECTORIES

    #print2['Xi'] = np.round(print2['POSITION_X_t2'].div(px), decimals=3)*(scale)  #PARTIAL TRAJECTORIES
    #print2['Yi'] = np.round(print2['POSITION_Y_t2'].div(px), decimals=3)*(scale)  #PARTIAL TRAJECTORIES

    kf22 = print2.loc[Tracks2,:]
    tracks_2 = print2.groupby('TRACK_ID2')

    img2 = Image.new("RGB", (L,L), "black")   # IMPORTANT to set dimensions of new image <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    draw2 = ImageDraw.Draw(img2)
    colors = itertools.cycle(('red','red'))
    #"Red", "Chartreuse", "DarkOrange", "DarkMagenta", "gold",'gray','navy','plum','deepskyblue','purple','maroon'
    for o, group in tracks_2:
        coords = tuple(zip(group.Xi, group.Yi))
        draw2.line(coords, width=25, fill=next(colors))  #to draw each trajectory with different colours
        draw2.ellipse([(coords[0]),(tuple([50+x for x in coords[0]]))], fill='yellow')
    #   draw.point(coords,fill='yellow')    #to draw each trajectory with the same points
    #   draw.line(coords,width=1, fill='red')    #to draw each trajectory with the same colour

    maxsize = (pk, pk)
    img2 = img2.resize(maxsize, Image.BICUBIC)
    #img2.show()
    traj2 = img2.save((file)+ "/TRACKS2.tif","tiff", transparency=0)



import os; os.system('say "Done!"')
    