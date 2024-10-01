from pywavesurfer import ws
import skimage as ski
import pandas as pd
import numpy as np
import scipy as sp
import stumpy
import os
from unityvr.preproc import logproc as lp
from unityvr.preproc import logproc as lp
from unityvr.viz import viz
from unityvr.analysis import posAnalysis, shapeAnalysis, fitting, utils


def importWaveSurfer(dirName, h5fileName, samplingFreq = 20e3):
    wsDict = ws.loadDataFile(filename=dirName+'/'+h5fileName, format_string='double' )
    wsDf = pd.DataFrame(wsDict[list(wsDict)[1]]['analogScans'].T, columns = list(wsDict[list(wsDict)[0]]['AIChannelNames'].astype('str')))
    
    #get wavesurfer time
    T = len(wsDf)/samplingFreq
    time = np.linspace(0,T,len(wsDf))
    wsDf['time'] = time
    return wsDf

'''def alignUnityWS(posDf, wsDf):
    #downsample wavesurfer file to unity frame rate
    wsDf_ds = wsDf.iloc[::int(np.round(np.nanmean(posDf['time'].diff())/np.nanmean(wsDf['time'].diff()),1))].copy().reset_index(drop=True)
    #pad downsampled wavesurfer with 0s and then use phase cross correlation to compute shift
    shift = ski.registration.phase_cross_correlation(np.hstack([wsDf_ds['Arena Heading'],
                                                        np.zeros(len(posDf)-len(wsDf_ds))]), (360-posDf['angle']).values)[0][0]
    if shift>0:
        raise ValueError('Wavesurfer file seems to have started before unity')
    else:
        posDf = posDf.iloc[-int(np.round(shift)):,:].copy().reset_index(drop=True)
    return posDf, wsDf_ds
    '''

def alignUnityWS(posDf, wsDf, clip=False):
    if len(wsDf)<len(posDf):
        distance_profile = stumpy.mass(wsDf['Arena Heading'], 360-posDf['angle'])
        shift = np.argmin(distance_profile)
        print(shift)
        if shift<0:
            raise ValueError('Wavesurfer file seems to have started before unity even though wavesurfer demarcates')
        else:
            posDf = posDf.iloc[-int(np.round(shift)):,:].copy().reset_index(drop=True)
    else:
        distance_profile = stumpy.mass(360-posDf['angle'], wsDf['Arena Heading'])
        shift = np.argmin(distance_profile)
        print("shift value: ", shift)
        if clip:
            if shift<0:
                raise ValueError('Unity file seems to have started before wavesurfer even though unity demarcates')
            else:
                wsDf_clipped = wsDf.iloc[int(np.round(shift)):int(np.round(shift))+len(posDf),:].copy().reset_index(drop=True)
                #wsDf_clipped.iloc[:len(posDf),:] #only the section of wavesurfer that is a match is returned
        else:
            wsDf_clipped = wsDf.copy()
    
    return posDf, wsDf_clipped, shift

def find_files_by_extension(dirName, extension1):
    file_list = []
    for root, dirs, files in os.walk(dirName):
        for file in files:
            if file.endswith(extension1):
                file_list.append(file)
                
    # Sort files by date created
    file_list.sort(key=lambda x: os.path.getctime(os.path.join(dirName, x)))
    return file_list

def mergeDfs(dirName, unityFiles, wsFiles, plot=False, kinematic_subject=False):
    if len(unityFiles)>len(wsFiles):
        fullDf = pd.DataFrame()
        fullVidDf = pd.DataFrame()
        wsDf = importWaveSurfer(dirName,wsFiles[0])
        cum_val_rot = 0;
        for i,f in enumerate(unityFiles):
            print(f)
            uvrTest = lp.constructUnityVRexperiment(dirName,f)
            #cum_val_rot = cum_val_rot + uvrTest.posDf['angle'].iloc[0]
            posDf = posAnalysis.position(uvrTest, 
                             rotate_by = None,
                             derive = True, #derive set to true adds 
                             #derived parameters like velocity and angle to the dataframe
                             plot=plot
                             #pass the following parameters to save the dataframe in a chosen directory
                             #,plotsave=False,saveDir=saveDir
                            )
            #account for weird example kinematic effect
            if kinematic_subject: 
                if len(posDf['time'].unique()) != len(posDf):
                    posDf = posDf.iloc[::2].reset_index(drop=True)
                posDf = posDf.fillna(0)
                
            vidDf = uvrTest.vidDf.copy()
            wsDf_ds = wsDf.iloc[::int(np.round(np.nanmedian(posDf['time'].diff(
            ))/np.nanmedian(wsDf['time'].diff()),1))].copy().reset_index(drop=True)
            
            posDf, wsDf_ds_clipped, shift = alignUnityWS(posDf, wsDf_ds, clip=True)
            posDf['laser'] = wsDf_ds_clipped['Laser Punishment'].values
            posDf['trial'] = i
            posDf['trial_name'] = f

            #ASSUMPTION: UNITY X AND Y START AT 0, if problem fixed on unity end, remove snippet from below
            #relies on sorting by time
            if 'x' in fullDf:
                posDf['x'] = posDf['x'] + fullDf['x'].iloc[-1]
                posDf['y'] = posDf['y'] + fullDf['y'].iloc[-1]
                posDf['frame'] = posDf['frame'] + fullDf['frame'].iloc[-1]

            #breaks independance
                vidDf['frame'] = vidDf['frame'] + fullDf['frame'].iloc[-1]
            
            fullDf = pd.concat([fullDf,posDf])
            fullVidDf = pd.concat([fullVidDf,vidDf])
        return fullDf.reset_index(drop=True), fullVidDf.reset_index(drop = True), wsDf_ds
        

def align_UVR_and_2p(diffs_idx_UVR, diffs_idx_2p, fpv=16):
    arenaDf = pd.DataFrame()
    arenaDf['Time Sec'] = diffs_idx_UVR
    arenaDf['Time Sec'] = arenaDf['Time Sec']/20000
    arenaDf['Unity Index'] = range(1,len(arenaDf)+1)

    two_pDf = pd.DataFrame()
    two_pDf['Time Sec'] = diffs_idx_2p
    two_pDf['Time Sec'] = two_pDf['Time Sec']/20000

    volDf = two_pDf.iloc[::fpv].reset_index(drop=True)
    volDf['Vol Index'] = range(1, len(volDf) +1)

    tempDf = pd.merge(arenaDf, volDf, how='outer', on='Time Sec').sort_values('Time Sec').reset_index(drop = True)
    tempDf['Unity Index'] = tempDf['Unity Index'].fillna(method='ffill', inplace=False)
    alignDf = tempDf.loc[~np.isnan(tempDf['Vol Index'])].reset_index(drop = True)

    return alignDf

def count_volumes(data, threshold=2.5, numFpv=16, clock = '2p'):
    crossings = 0
    previous_value = data[0]
    diffs_index = []
    i = 0

    for current_value in data[1:]:
        i += 1
        if clock != '2p': 
            if previous_value >= threshold and current_value < threshold or previous_value <= threshold and current_value > threshold:
                crossings += 1
                diffs_index.append(int(i-1))  # Append instead of using index assignment
            previous_value = current_value
        else:
            if previous_value >= threshold and current_value < threshold:
                crossings += 1
                diffs_index.append(int(i-1))  # Append instead of using index assignment
            previous_value = current_value

    return ['Frames: ' + str(crossings) + '  Volumes: ' + str(crossings / numFpv), np.diff(diffs_index), [int(i) for i in diffs_index]]

    
    