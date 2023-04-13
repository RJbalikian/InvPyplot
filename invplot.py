import csv
import pathlib 
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ipywidgets as widgets
import scipy.interpolate
import matplotlib


def ingest_inv(inv_file, verbose=True):
    """Function to ingest inversion file and get key points (row numbers) in the file

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        The res2Dinv .inv file to work with.

    Returns
    -------
    inv_dict : dict
        Dictionary containing the important locations in the file
    """
    if isinstance(inv_file, pathlib.PurePath):
        pass
    else:
        inv_file = pathlib.Path(inv_file)
    
    fileHeader = []
    iterationStartRowList = []
    layerRowList = []
    layerDepths = []
    noLayerRow = -1
    blockRow = -1
    layerRow = -1
    layerInfoRow = -1
    resistDF = pd.DataFrame()
    dataList = []
    noPoints = []
    calcResistivityRowList = []
    refResistRow=-1
    topoDataRow = -1
    iterationsInfoRow = -1

    with open(str(inv_file)) as datafile: 
        filereader = csv.reader(datafile)
        for row in enumerate(filereader):
            startLayer = 0
            endLayer = 0
            lay = -1
            #print(row[0])
            if row[0] <= 8:
               
                if len(row[1])>1:
                    fileHeader.append(row[1][0]+', '+row[1][1])
                    continue
                else:
                    fileHeader.append(row[1][0].strip())
                    continue
                
            if 'NUMBER OF LAYERS' in str(row[1]):
                noLayerRow = row[0]+1
                continue
            if row[0] == noLayerRow:
                noLayers = int(row[1][0])
                layerList = np.linspace(1,noLayers, noLayers)
                continue
            
            if 'NUMBER OF BLOCKS' in str(row[1]):
                blockRow = row[0]+1
                continue
            if row[0]==blockRow:
                noBlocks = int(row[1][0])
                continue

            if 'ITERATION' in str(row[1]):
                iterationStartRowList.append(row[0]) #Add row of iteration to iterationStartRowList
                continue

            if 'LAYER ' in str(row[1]):
                iterInd = len(iterationStartRowList)-1
                if iterInd > len(layerRowList)-1:
                    layerRowList.append([row[0]])
                else:
                    layerRowList[iterInd].append(row[0])
                layerInfoRow = row[0]+1
                continue
            if row[0]==layerInfoRow:
                noPoints.append(int(row[1][0].strip()))
                layerDepths.append(row[1][1].strip())
                continue
            
            if 'CALCULATED APPARENT RESISTIVITY' in str(row[1]):
                calcResistivityRowList.append(row[0])
                continue
            
            if 'Reference resistivity is' in str(row[1]):
                refResistRow = row[0]+1 
                continue
            
            if row[0]==refResistRow:
                refResist = float(row[1][0].strip())
                continue
            
            if 'TOPOGRAPHICAL DATA' in str(row[1]):
                topoDataRow = row[0]
                continue
            if row[0]==topoDataRow+2:
                noTopoPts = int(row[1][0].strip())
                continue

            if 'COORDINATES FOR ELECTRODES' in str(row[1]):
                electrodeCoordsRow = row[0]
                continue

            if 'Shift matrix' in str(row[1]):
                shiftMatrixRow = row[0]
                continue

            if 'Blocks sensitivity and uncertainity values (with smoothness constrain)' in str(row[1]):
                sensAndUncertainRow = row[0]
                continue

            if 'Error Distribution' in str(row[1]):
                errorDistRow = row[0] #no of data points
                continue

            if 'Total Time' in str(row[1]):
                iterationsInfoRow=row[0]
                iterDataList = []
                continue
            
            if iterationsInfoRow > 1:
                if row[1] == []:
                    print('   ')
                    noIterations = row[0]-iterationsInfoRow-1
                    break
                iterDataList.append(row[1][0].split())

    layerDepths = layerDepths[0:noLayers]
    layerDepths[noLayers-1] = float(layerDepths[(noLayers-2)])+(float(layerDepths[noLayers-2])-float(layerDepths[noLayers-3]))
    layerDepths = [float(x) for x in layerDepths]

    noPoints = noPoints[0:noLayers]

    keyList=['Name', 'NomElectrodeSpacing', 'ArrayCode', 'ProtocolCode', 'MeasureHeader', 'MeasureType', 'NoDataPoints','DistanceType','FinalFlag']

    global fileHeaderDict
    fileHeaderDict = dict(zip(keyList, fileHeader))
    noDataPoints = int(fileHeaderDict['NoDataPoints'])
    iterationDF = pd.DataFrame(iterDataList, columns=['Iteration',  'Time for this iteration', 'Total Time', '%AbsError'])
    iterationDF = iterationDF.apply(pd.to_numeric)

    if verbose:
        print(iterationDF)
        fig1, ax1 = plt.subplots(1)
        iterationDF.plot('Iteration','%AbsError',figsize=(2,2), ax=ax1, c='k')
        iterationDF.plot('Iteration','%AbsError',figsize=(2,2),kind='scatter', ax=ax1, c='k')
        ax1.set_title(inv_file.stem)
        ax1.set_xticks(np.arange(0,iterationDF['Iteration'].max()+1))
        ax1.get_legend().remove()
        plt.show(fig1)

    inv_dict = {
     'fileHeader':fileHeader,
     'iterationStartRowList':iterationStartRowList,
     'layerRowList':layerRowList, 
     'layerDepths':layerDepths,
     'noLayerRow':noLayerRow,
     'blockRow':blockRow ,
     'layerRow':layerRow ,
     'layerInfoRow':layerInfoRow ,
     'resistDF':resistDF,
     'dataList':dataList,
     'noPoints':noPoints,
     'calcResistivityRowList':calcResistivityRowList ,
     'refResistRow':refResistRow,
     'topoDataRow':topoDataRow,
     'iterationsInfoRow':iterationsInfoRow,
     'iterationDF':iterationDF,
     'noIterations':noIterations,
     'noDataPoints':noDataPoints,
     'shiftMatrixRow':shiftMatrixRow,
     'electrodeCoordsRow':electrodeCoordsRow,
     'noTopoPts':noTopoPts,
     'sensAndUncertainRow':sensAndUncertainRow,
     'noModelBlocks':[],
     'errorDistRow':errorDistRow,
     'fileHeaderDict':fileHeaderDict}

    return inv_dict
