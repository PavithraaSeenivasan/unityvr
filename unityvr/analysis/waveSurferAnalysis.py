from pywavesurfer import ws
import skimage as ski
import pandas as pd
import numpy as np
import scipy as sp


def importWaveSurfer(dirName, h5fileName, samplingFreq = 20e3):
    wsDict = ws.loadDataFile(filename=dirName+'/'+h5fileName, format_string='double' )
    wsDf = pd.DataFrame(wsDict[list(wsDict)[1]]['analogScans'].T, columns = list(wsDict[list(wsDict)[0]]['AIChannelNames'].astype('str')))
    
    #get wavesurfer time
    T = len(wsDf)/samplingFreq
    time = np.linspace(0,T,len(wsDf))
    wsDf['time'] = time
    return wsDf

def alignPosDf(posDf, wsDf):
    #downsample wavesurfer file to unity frame rate
    wsDf_ds = wsDf.iloc[::int(np.round(np.nanmean(posDf['time'].diff())/np.nanmean(wsDf['time'].diff()),1))].copy().reset_index(drop=True)
    #pad downsampled wavesurfer with 0s and then use phase cross correlation to compute shift
    shift = ski.registration.phase_cross_correlation(np.hstack([wsDf_ds['Arena Heading'],
                                                        np.zeros(len(posDf)-len(wsDf_ds))]), (360-posDf['angle']).values)[0][0]
    if shift>0:
        raise ValueError('Wavesurfer file seems to have started before unity')
    else:
        posDf = posDf.iloc[shift:].copy().reset_index(drop=True)
    return posDf
        

    
    