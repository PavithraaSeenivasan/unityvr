import xarray as xr

from os.path import sep, exists
from os import mkdir, makedirs, getcwd

from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter
import skimage as ski

import numpy as np
import pandas as pd
import json
from dataclasses import dataclass, asdict
import ast

from matplotlib import pyplot as plt

def deriveTexVals(texDf, 
              std_filter = 3, #3*std deviation filter for removing large jumps
              diskSize = 5, #morphological disk
              round=-1, #rounding to the nearest 10th in deg/sec
              screenAboveFly = 32, #for pentagonal display with each screen dimension = 9.5*5.8 cm, in degs, +ve
              screenBelowFly = 60 #in degs, +ve
             ):
    texDf = texDf.copy()
    texDf['stimAngle'] = (-texDf['azimuth'].values)%360-180 #convert to -180 to 180 left handed convention
    vel = texDf['stimAngle'].diff().values #left handed convention
    vel[np.abs(vel)>(np.nanmean(vel)+std_filter*np.nanstd(vel))] = 0 #remove large jumps
    texDf['stimVel'] = np.round(ski.morphology.closing(ski.morphology.opening(vel,np.ones(diskSize)),
                                                        np.ones(diskSize))/np.nanmedian(texDf['time'].diff()),round)
    
    #apply morphological operation to remove unity noise and round off to the nearest 10th
    A = 1/(np.tan(screenAboveFly*np.pi/180)-np.tan(-screenBelowFly*np.pi/180))
    B = -np.tan(-screenBelowFly*np.pi/180)
    elevationToDegs = lambda e : np.arctan(e/A-B)*180/np.pi
    texDf['elevationDegs'] = elevationToDegs(texDf['elevation'].values).round(0)
    texDf['stimSpeed'] = np.abs(texDf['stimVel'])
    texDf['stimDir'] = np.sign(texDf['stimVel'])
    texDf['behindScreen'] = np.abs(texDf['stimAngle'])>=(180-36);
    return texDf

def convertTextureVals(texDf, RF=True):
    if RF: 
        #elevation was mapped
        texDf['elevation'] = np.round(1-(texDf.ytex % 1),2)
    texDf.xtex = texDf.xtex - texDf.xtex[0]
    xtexpos = texDf.xtex.values.copy()
    xtexpos[xtexpos<0] = 1+xtexpos[xtexpos<0]
    texDf['azimuth'] = xtexpos*360
    texDf['sweepdir'] = np.sign(texDf.xtex) #right handed convention
    return texDf

def deriveVidVals(uvrDat, movieFolderPath, imageFile = 'stimGenDf.csv', sceneFile='scene1DArray.npy', shift = -90):
    movieFolder = uvrDat.vidDf['img'].str.split(r'\\').str.get(-2).unique()[-1]
    moviePath = movieFolderPath + movieFolder
    uvrDat.vidDf['filename'] = uvrDat.vidDf['img'].str.split(r'\\').str.get(-1)
    if imageFile is not None:
        stimGenDf = pd.read_csv(sep.join([moviePath,imageFile]),index_col=0)
        uvrDat.vidDf = pd.merge(uvrDat.vidDf, stimGenDf, on=['filename'])
    columnsToKeepVid = list(uvrDat.vidDf.columns)
    uvrDat.vidDf = pd.merge(uvrDat.posDf,uvrDat.vidDf,on = ['frame'],how='left').drop(columns='time_y').rename(columns={'time_x':'time'}).ffill()[columnsToKeepVid]
    if sceneFile is not None:
        uvrDat.sceneArray = np.load(sep.join([moviePath,sceneFile]))
        uvrDat.sceneArray = np.roll(uvrDat.sceneArray[:,:], shift=int(np.round(uvrDat.sceneArray.shape[-1]/360*shift)))
    return uvrDat