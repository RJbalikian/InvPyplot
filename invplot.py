import csv
import pathlib 
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate
import matplotlib

def autoplot(inv_file, iteration, return_dict=False, **kwargs):
    """Function to run all intermedaite functions and resinv_plot to simply read and plot everything in one call.

    Parameters
    ----------
    inv_file : str or pathlib.PurePath object
        Filepath to .inv file of interest. The .inv file should be one generated from Res2DInv.
    iteration : int or list or either str {':', 'all'}
        Integer or list of integers indicating which iteration of the .inv result to use for plotting. If list, all will be plotted separately. If ':' or 'all', will plot all iterations successively.
    return_dict : bool, optional
        Whether to return results as a dictionary, by default False
    **kwargs
        Other keyword arguments may be read into autoplot. These are read in as **kwargs to either resinv_plot() or matplotlib.pyplot.imshow via the resinv_plot function. See documentation for resinv_plot for available parameters for resinv_plot.

    Returns
    -------
    inv_dict : dict
        Dictionary containing all input parameters and data generated along the way, including the output figures and axes

    """
    if isinstance(inv_file, pathlib.PurePath):
        pass
    else:
        inv_file = pathlib.Path(inv_file)

    inv_dict = ingest_inv(inv_file, verbose=False, show_iterations=False)
    inv_dict = read_inv_data(inv_file=inv_file, inv_dict=inv_dict)

    allIterList = [':', 'all']
    if type(iteration) is int:
        iteration = [iteration]
    elif iteration.lower() in allIterList:
        iteration = inv_dict['iterationDF'].Iteration.tolist()

    resinv_params_list = ['inv_dict', 'colMap', 'cBarFormat', 'cBarLabel', 'cBarOrientation', 'cMin', 'cMax', 
                          'griddedFt', 'griddedM', 'title', 'normType', 'primaryUnit', 'showPoints','whichTicks', 
                          'figsize', 'dpi', 'reverse', 'tight_layout', 'savefig', 'saveformat']
    resinv_kwargs = {}
    imshow_kwargs = {}
    for key, value in kwargs.items():
        if key in resinv_params_list:
            resinv_kwargs[key] = value
        else:
            imshow_kwargs[key] = value

    iterIndList = []
    iterNoList = []
    figList = []
    axList = []

    for i in iteration:
        iterNo = i
        inv_dict['iterationNo'] = iterNo
        iterInd = inv_dict['iterationDF'][inv_dict['iterationDF'].Iteration==i].index.tolist()[0]
        inv_dict['iterationInd'] = iterInd
        inv_dict = read_inv_data_other(inv_file=inv_file, inv_dict=inv_dict, iteration_no=iterNo)
        inv_dict = read_error_data(inv_file=inv_file, inv_dict=inv_dict)
        inv_dict = get_resistivitiy_model(inv_file=inv_file, inv_dict=inv_dict)
        fig, ax = resinv_plot(inv_dict=inv_dict, imshow_kwargs=imshow_kwargs, **kwargs)

        iterIndList.append(i)
        iterNoList.append(inv_dict['iterationDF'].loc[iterInd, 'Iteration'])
        figList.append(fig)
        axList.append(ax)

    inv_dict['iterationNo'] = iterIndList
    inv_dict['iterationInd'] = iterNoList
    inv_dict['fig'] = figList
    inv_dict['ax'] = axList
    
    if return_dict:
        return inv_dict
    return

#Function that performs all the actual plotting
def resinv_plot(inv_dict, colMap='nipy_spectral', cBarFormat ='%3.0f', cBarLabel='Resistivity (ohm-m)', cBarOrientation='horizontal', cMin=None, cMax=None, griddedFt=[False,False], griddedM=[False,False], title=None, normType='log', primaryUnit='m', showPoints=False,whichTicks='major', figsize=None, dpi=None, reverse=False, tight_layout=True, savefig=False, saveformat='png', imshow_kwargs=None, **kwargs):
    """Function to pull everything together and plot it nicely.

    It is recommended to use the autoplot function rather than resinv_plot directly, since using autoplot() incorporates all the setup needed to create the input dictionary keys/values correctly.

    Parameters
    ----------
    inv_dict : dict
        Dictionary of inversion results generated from previous steps
    colMap : str, optional
        Colormap, any acceptable from matplotlib, by default 'nipy_spectral'
    cBarFormat : str, optional
        Format string for colorbar tick labels, by default '%3.0f'
    cBarLabel : str, optional
        Colorbar label, by default 'Resistivity (ohm-m)'
    cBarOrientation : str {'horizonta', 'vertical'}, optional
        Orientation of the colorbar, by default 'horizontal'
    cMin : float, optional
        Minimum of colorbar/colormap, by default None, which uses the minimum value of the dataset.
    cMax : float, optional
        Maximum of colorbar/colormap, by default None, which uses the maximum value of the dataset.
    griddedFt : list, optional
        Whether to show gridlines on the feet ticks, first position is x, second position is y, by default [False,False]
    griddedM : list, optional
        Whether to show gridlines on the meter tickes, first position is x, second posistion is y, by default [False,False]
    title : str, optional
        String to show as the title, if desired to set manually, by default None, which shows the filename as the title
    normType : str {'log', 'linear'}, optional
        Normalization type, by default 'log'. Determines whether matplotlib.colors.LogNorm or matplotlib.colors.Normalize is used for colormap.
    primaryUnit : str {'m', 'ft'}, optional
        Whether to display meters or feet as primary unit (this determines which unit is larger on the axis and is on the left and top), by default 'm'
    showPoints : bool, optional
        Whether to show the datapoints used for interpolation, by default False
    whichTicks : str {'major', 'minor', 'both'}, optional
        If griddedFt or griddedM has any True, this determines whether major, minor, or both gridlines are used; by default 'major'.
    figsize : tuple, optional
        Tuple (width, height) of the figsize, read into plt.rcParams['figure.figsize'], by default None.
    dpi : int or float, optional
        Resolution (dots per square inch) of final figure, read into plt.rcParams['figure.dpi'], by default None.
    reverse : bool, optional
        Whether to display the data in reverse (flipped along x) of what is read into from .inv file, by default False
    tight_layout : bool, optional
        If true, calls fig.tight_layout(). Otherwise, tries to maximize space on the figure using plt.subplots_adjust, by default True
    savefig : bool, optional
        If False, will not save figure. Otherwise, calls plt.savefig() and the value of this parameter will be used as the output filepath, by default False.
    saveformat : str, optional
        Read into plt.savefig(format) paramater, by default 'png'.

    Returns
    -------
    dict
        Returns existing inv_dict input, but with added keys of ['fig'] and ['ax'] containing a list of the fig and ax objects (list, since multiple iterations can be done at once)
    """
    if title is None:
        title = inv_dict['inv_file_Path'].stem
    
    if 'figure.dpi' not in list(inv_dict.keys()):
        inv_dict['figure.dpi'] = 250
    if 'figure.figsize' not in list(inv_dict.keys()):
        inv_dict['figure.figsize'] = (12,5)

    x = inv_dict['resistModelDF']['x'].copy()
    z = inv_dict['resistModelDF']['zElev'].copy()
    v = inv_dict['resistModelDF']['Data'].copy()

    if figsize is None:
        plt.rcParams['figure.figsize'] = inv_dict['figure.figsize']
    else:
        inv_dict['figure.figsize'] = figsize
        plt.rcParams['figure.figsize'] = figsize

    if dpi is None:
        plt.rcParams['figure.dpi'] = inv_dict['figure.dpi']
    else:
        inv_dict['figure.dpi'] = dpi
        plt.rcParams['figure.dpi'] = dpi
    
    maxXDist = max(np.float_(inv_dict['electrodes']))

    if cMin is None:
        cMin = inv_dict['resistModelDF']['Data'].min()
    if cMax is None:
        cMax = inv_dict['resistModelDF']['Data'].max()
    

    for i in enumerate(v):
        v[i[0]] = abs(float(i[1]))

    xi, zi = np.linspace(min(x), max(x), int(max(x))), np.linspace(min(z), max(z), int(max(z)))
    xi, zi = np.meshgrid(xi, zi)

    vi = scipy.interpolate.griddata((x, z), v, (xi, zi))#, method='linear')

    ptSize = round(100 / maxXDist * 35, 1)

    fig, axes = plt.subplots(1)
    cmap = matplotlib.cm.binary
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0,1,cmap.N)
    my_cmap = matplotlib.colors.ListedColormap(my_cmap)

    vmax98 = np.percentile(v, 98)
    vmin2 = np.percentile(v, 2)
    minx = min(x)
    maxx = max(x)
    minz = min(z)
    maxz = max(z)

    vmax = cMax
    vmin = cMin

    #if cMax >= resistModelDF['Data'].max():
    #  vmax = vmax98
    #else:
    #  vmax = cMax
    #if cMin <= resistModelDF['Data'].min():
    #  vmin = vmin2
    #else:
    #  vmin = cMin
    #cbarTicks = np.arange(np.round(vmin,-1),np.round(vmax-1)+1,10) 

    arStep = np.round((vmax-vmin)/10,-1)
    cbarTicks = np.arange(np.round(vmin, -1), np.ceil(vmax/10)*10,arStep)

    #Get default values or kwargs, depending on if kwargs have been used
    if 'norm' in imshow_kwargs.keys():
        norm = imshow_kwargs['norm']
    else:
        if normType=='log':
            if vmin <= 0:
                vmin = 0.1
            norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax)
            #cBarFormat = '%.1e'
            #cbarTicks = np.logspace(np.log10(vmin),np.log10(vmax),num=10)
        else:
            norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)

    #im = self.axes.imshow(vi, vmin=vmin, vmax=vmax, origin='lower',
    if 'extent' in imshow_kwargs.keys():
        extent = imshow_kwargs['extent']
        imshow_kwargs.pop('extent', None)
    else:
        extent = [minx, maxx, minz, maxz]
        
    if 'aspect' in imshow_kwargs.keys():
        aspect = imshow_kwargs['aspect']
        imshow_kwargs.pop('aspect', None)
    else:
        aspect = 'auto'
        
    if 'cmap' in imshow_kwargs.keys():
        cmap = imshow_kwargs['cmap']
        imshow_kwargs.pop('cmap', None)
    else:
        cmap=colMap   
        
    if 'interpolation' in imshow_kwargs.keys():
        interp = imshow_kwargs['interpolation']
        imshow_kwargs.pop('interpolation', None)
    else:
        interp='spline36'   
    
    im = axes.imshow(vi, origin='lower',
                extent=extent,
                aspect=aspect,
                cmap =cmap,
                norm = norm,
                interpolation=interp, **imshow_kwargs)
    f, a = __plot_pretty(inv_dict, x,z,v,fig=fig,im=im,ax=axes,colMap=colMap,cMin=cMin,cMax=cMax, 
                       gridM=griddedM, gridFt=griddedFt, primaryUnit=primaryUnit, t=title, tight_layout=tight_layout,
                       cbarTicks=cbarTicks,cBarFormat=cBarFormat,cBarLabel=cBarLabel,cBarOrient=cBarOrientation,
                       showPoints=showPoints, norm=norm, whichTicks=whichTicks, reverse=reverse)

    plt.show(fig)
    if savefig is not False:
        plt.savefig(savefig, format=saveformat, facecolor='white')

    plt.close(f)
    return f, a

#Function to ingest inv file and find key information for later
def ingest_inv(inv_file, verbose=True, show_iterations=True):
    """Function to ingest inversion file and get key points (row numbers) in the file

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        The res2Dinv .inv file to work with.
    verbose : bool, default=True
        Whether to print results to terminal. Here, prints a pandas dataframe with information about iterations.
    show_iterations : bool, default=True
        Whether to show a matplotlib plot with the iteration on x axis and percent error on y axis.

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
    if show_iterations:
        fig1, ax1 = plt.subplots(1)
        iterationDF.plot('Iteration','%AbsError',figsize=(3,3), ax=ax1, c='k')
        iterationDF.plot('Iteration','%AbsError',figsize=(3,3),kind='scatter', ax=ax1, c='k')
        ax1.set_title(inv_file.stem)
        ax1.set_xticks(np.arange(0,iterationDF['Iteration'].max()+1))
        ax1.get_legend().remove()
        plt.show(fig1)

    inv_dict = {
        'inv_file_Path':inv_file,
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

#Input Data
def read_inv_data(inv_file, inv_dict, startRow=9):
    """Function to do initial read of .inv file. 
    
    This data does not change with iteration, as in later function. This function should be readafter ingest_inv, using the output from that as inv_dict.
     
    

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        Filepath of .inv file
    inv_dict : dict
        Dictionary (output from ingest_inv) containing information about where data is located in the .inv file
    startRow : int, optional
        Where to start the read. This is, by default 9, which works well for .inv files tested.

    Returns
    -------
    _type_
        _description_
    """
    noDataPoints = inv_dict['noDataPoints']
    if isinstance(inv_file, pathlib.PurePath):
        pass
    else:
        inv_file = pathlib.Path(inv_file)

    import csv
    with open(inv_file) as datafile: 
        filereader = csv.reader(datafile)
        start = 0
        inDataList = []
        for row in enumerate(filereader):
            if row[0] < startRow:
                continue
            elif row[0] < startRow+noDataPoints:
                inDataList.append(re.sub('\s+',' ',row[1][0]).split(' '))
            else:
                break
    inDF = pd.DataFrame(inDataList)
    if startRow == 9:
        inDF.drop([0],inplace=True,axis=1)
    inDF.astype(np.float64)
    inDF.columns=['NoElectrodes', 'A(x)', 'A(z)', 'B(x)', 'B(z)', 'M(x)', 'M(z)', 'N(x)', 'N(z)', 'Data']
    inv_dict['resistDF']  = inDF
    return inv_dict

#Read other important inversion data
def read_inv_data_other(inv_file, inv_dict, iteration_no=None):
    """Function to read inversion data. 

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        Filepath to .inv file of interest.
    inv_dict : dict
        Dictionary contianing outputs from ingest_inv and read_in_data.
    iteration_no : int
        Iteration number of interest.

    Returns
    -------
    dict
        Dictionary with more information appended to input dictionary
    """
    if iteration_no is None:
        print('Please run read_inv_data_other again and specify iteration by setting iteration_no parameter equal to integer.')
        return
    #Extract needed variables from dict
    iterationInd = inv_dict['iterationDF'][inv_dict['iterationDF'].Iteration==iteration_no].index.tolist()[0]
    inv_dict['iterationNo'] = iteration_no
    inv_dict['iterationInd'] = iterationInd


    invDF = inv_dict['resistDF']
    layerRowList = inv_dict['layerRowList']
    noPoints = inv_dict['noPoints']
    layerDepths = inv_dict['layerDepths']
    shiftMatrixRow = inv_dict['shiftMatrixRow']
    electrodeCoordsRow = inv_dict['electrodeCoordsRow']
    topoDataRow = inv_dict['topoDataRow']
    noTopoPts = inv_dict['noTopoPts']
    sensAndUncertainRow = inv_dict['sensAndUncertainRow']

    #Get Electrodes
    electrodes= pd.concat([invDF['A(x)'],invDF['B(x)'], invDF['M(x)'], invDF['B(x)'], invDF['N(x)']],ignore_index=True)
    electrodes.reset_index(inplace=True, drop=True)
    electrodes = electrodes.unique()
    inv_dict['electrodes'] = electrodes

    #ElectrodeCoordinates
    noModelElects = shiftMatrixRow-electrodeCoordsRow-1
    electrodeCoordsDF = pd.read_table(inv_file,skiprows=electrodeCoordsRow,nrows=noModelElects, sep='\s+')
    electrodeCoordsDF.dropna(axis=1,inplace=True)
    electrodeCoordsDF.columns=['xDist','RelElevation']
    electrodeCoordsDF['ElectrodeNo'] = electrodeCoordsDF.index+1
    inv_dict['electrodeCoordsDF'] = electrodeCoordsDF

    #Topographical Data
    topoDF = pd.read_table(inv_file,skiprows=topoDataRow+2,nrows=noTopoPts, sep='\s+')
    topoDF.reset_index(inplace=True)
    topoDF.columns=['xDist','Elevation']
    topoDF['ElectrodeNo'] = topoDF.index+1
    inv_dict['topoDF'] = topoDF

    #Resistivity Model
    resistModelDF = pd.DataFrame()
    for r in enumerate(layerRowList[iterationInd]):
        layerDepth = layerDepths[r[0]]
        noPtsInLyr = noPoints[r[0]]
        currDF = pd.read_table(inv_file,skiprows=r[1]+1, nrows=noPtsInLyr,sep=',')
        currDF.columns=['ElectrodeNo','Data']
        currDF['z'] = layerDepth
        resistModelDF= pd.concat([resistModelDF,currDF],ignore_index=True).copy()
    resistModelDF.reset_index(inplace=True, drop=True)
    noModelBlocks=resistModelDF.shape[0]
    inv_dict['resistModelDF'] = resistModelDF
    inv_dict['noModelBlocks'] = noModelBlocks

    #Shift Matrix
    shiftMatrixDF = pd.read_table(inv_file,skiprows=shiftMatrixRow+1,nrows=noModelElects, sep='\s+',header=None,index_col=0)
    shiftMatrixDF.dropna(axis=1,inplace=True)
    for c in shiftMatrixDF:
        shiftMatrixDF.rename(columns={c:'Layer'+str(int(c)-1)},inplace=True)
    inv_dict['shiftMatrixDF'] = shiftMatrixDF

    #Sensitivity
    sensDF = pd.read_table(inv_file,skiprows=sensAndUncertainRow+3,nrows=noModelBlocks, sep='\s+',header=None,index_col=0)
    sensDF.dropna(axis=1,inplace=True)
    sensDF.reset_index(inplace=True)
    sensDF.columns=['BlockNo','Sensitivity', '%ApproxUncertainty']
    inv_dict['sensDF'] = sensDF

    return inv_dict

#Error Distribution
def read_error_data(inv_file, inv_dict):
    """Function to read data pertaining to model error

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        Filepath to .inv file of interest
    inv_dict : dict 
        Dictionary containing cumulative output from ingest_inv, read_inv_data, and read_inv_data_other

    Returns
    -------
    dict
        Ouptput dictionary containing all information from read_error_data and previous functions.
    """
    import csv
    noDataPoints = inv_dict['noDataPoints']
    startRow = inv_dict['errorDistRow']+1

    with open(inv_file) as datafile: 
        filereader = csv.reader(datafile)
        inDataList = []
        for row in enumerate(filereader):
            if row[0] < startRow:
                continue
            elif row[0] < startRow+noDataPoints:
                newrow = row[1]
                newrow.append(newrow[4][8:])
                newrow[4] = newrow[4][0:8]
                newrow = [x.strip() for x in newrow]
                inDataList.append(newrow)
            else:
                break
    inDF = pd.DataFrame(inDataList)
    inv_dict['errDistDF'] = inDF

    errDistDF = inv_dict['errDistDF']

    colList=['xDist?','nFactor?','Measure1','Measure2','PercentError','MoreStacks','AvgMeasure']
    for i in range(0,5):
        errDistDF[i]=errDistDF[i].astype(np.float64)
        errDistDF.rename(columns={i:colList[i]},inplace=True)
    errDistDF['AvgMeasure'] = (errDistDF['Measure1']+errDistDF['Measure2'])/2
    inv_dict['errDistDF'] = errDistDF
    return inv_dict

#Interpolate between two points, simple
def map_diff(xIn, x1,x2,y1,y2):
    """Simple, linear interpolation between two points

    Parameters
    ----------
    xIn : float, int, or numeric
        X Location at which y-value is desired. This should fall between (or be equal to) x1 and x2
    x1 : float, int, or numeric
        Initial X location, with known y-value y1 
    x2 : float, int, or numeric
        Second X location, with known y-value y2
    y1 : float, int, or numeric
        Known y-value at x1
    y2 : float, int, or numeric
        Known y-value at x2

    Returns
    -------
    float
        Y-value that is proportionally scaled based on the xIn relative distance to x1 and x2
    """
    if x1==xIn:
        yOut=y1
    elif x2==xIn:
        yOut = y2
    else:
        totXDiff = x2-x1
        percXDiff = (xIn-x1)/totXDiff
        totYDiff = y2-y1
        yOut = y1 + totYDiff*percXDiff
    return yOut

#Function to get resistivity model as pandas dataframe
def get_resistivitiy_model(inv_file, inv_dict):
    """Function to read the resistivity model, given the iteration of interest.

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        Filepath to .inv file of interest
    inv_dict : dict
        Dictionary containing cumulative output from invest_inv, read_inv_data, read_inv_data_other, and read_error_data.

    Returns
    -------
    dict
        Dictionary with resistivity model contained in new key:value pair of inv_dict
    """
    resistModelDF = pd.DataFrame()
    layerRowList= inv_dict['layerRowList']
    iterationInd= inv_dict['iterationInd']
    layerDepths= inv_dict['layerDepths']
    noPoints= inv_dict['noPoints']
    electrodeCoordsDF= inv_dict['electrodeCoordsDF']
    topoDF= inv_dict['topoDF']

    for r in enumerate(layerRowList[iterationInd]):
        layerDepth = layerDepths[r[0]]
        noPtsInLyr = noPoints[r[0]]
        currDF = pd.read_table(inv_file, skiprows=r[1]+1, nrows=noPtsInLyr,sep=',')
        currDF.iloc[currDF.shape[0]-1,1] =  currDF.iloc[currDF.shape[0]-1,0]
        currDF.iloc[currDF.shape[0]-1,0] = currDF.iloc[currDF.shape[0]-2,0]+1
        currDF.columns=['ElectrodeNo','Data']
        currDF['zDepth'] = layerDepth

        for i in currDF.index:
            lowerElecNo = currDF.loc[i,'ElectrodeNo']#-1
            elecInd = electrodeCoordsDF.loc[electrodeCoordsDF['ElectrodeNo']==lowerElecNo].index.values[0]
            currDF.loc[i,'x'] = (electrodeCoordsDF.loc[elecInd,'xDist'] + electrodeCoordsDF.loc[elecInd+1,'xDist'])/2
            for xT in enumerate(topoDF['xDist']):
                if xT[1] < currDF.loc[i,'x']:
                    continue
                else:
                    topoX1 = topoDF.loc[xT[0]-1,'xDist']
                    topoX2 = topoDF.loc[xT[0],'xDist']
                    topoZ1 = topoDF.loc[xT[0]-1,'Elevation']
                    topoZ2 = topoDF.loc[xT[0],'Elevation']
                    break
            currDF.loc[i,'zElev'] = map_diff(currDF.loc[i,'x'],topoX1, topoX2, topoZ1, topoZ2)-currDF.loc[i,'zDepth']
        if r[0] == 0:
            surfDF = currDF.copy()
            surfDF['zElev'] = surfDF.loc[:,'zElev']+surfDF.loc[:,'zDepth']
            surfDF['zDepth'] = 0
            resistModelDF = pd.concat([resistModelDF, surfDF], ignore_index=True)
            resistModelDF.reset_index(inplace=True, drop=True)
        resistModelDF = pd.concat([resistModelDF, currDF], ignore_index=True)
        resistModelDF.reset_index(inplace=True, drop=True)
    inv_dict['resistModelDF'] = resistModelDF
    return inv_dict

#Helper function for __plot_pretty
def __label_plot(fig, ax, gridM, gridFt, whichTicks, pUnit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims, yLims, t):
    """See __plot_pretty, as all input parameters are derived from there"""
    #matplotlib.rc('font', family='sans-serif')   
    #matplotlib.rc('font', serif='Helvetica') 
    #matplotlib.rc('text', usetex='false')   
    #fontName = {'fontname':'Helvetica'}
    
    plt.title(t)
    ax.set_title(t, fontsize=20)
    ax.set_xlabel('Distance ['+sUnit+']', fontsize = 12)
    ax.set_ylabel('Elevation ['+pUnit+']', fontsize = 14)
    ax.set_xticks(sUnitXLocs)
    ax.set_yticks(pUnitYLocs)
    ax.set_xticklabels(sUnitXLabels,fontsize=12)
    ax.set_yticklabels(pUnitYLabels,fontsize=12)
    ax.set_yticks(pUnitYLocs)
    ax.set_ylim(yLims)
    ax.set_xlim(xLims)
    ax.minorticks_on()

    ax2=ax.twiny()
    ax2.set_xticks(pUnitXLocs)
    ax2.set_xticklabels(pUnitXLabels,fontdict={'fontsize':14})
    ax2.set_xlabel('Distance ['+pUnit+']',fontsize=14)
    ax2.minorticks_on()
    ax2.set_xlim(xLims)

    ax3=ax2.twinx()
    ax3.set_yticks(sUnitYLocs)
    ax3.set_yticklabels(sUnitYLabels,fontdict={'fontsize':10})
    ax3.set_ylim(yLims)
    ax3.set_ylabel('Elevation ['+sUnit+']',fontsize=14, rotation=270, labelpad = 20)
    ax3.minorticks_on()
    

    if pUnit=='m':
        if gridM[0]:
           ax2.grid(axis='x',alpha=0.5, c='k', which=whichTicks)
        if gridM[1]:
            ax.grid(axis='y',alpha=0.5, c='k', which=whichTicks)
        if gridFt[0]:
            ax.grid(axis='x',alpha=0.5, c='k', which=whichTicks)
        if gridFt[1]:
            ax3.grid(axis='y',alpha=0.5, c='k', which=whichTicks)
    else:
        if gridFt[0]:
            ax2.grid(axis='x',alpha=0.5, c='k', which=whichTicks)
        if gridFt[1]:
            ax.grid(axis='y',alpha=0.5, c='k', which=whichTicks)
        if gridM[0]:
            ax.grid(axis='x',alpha=0.5, c='k', which=whichTicks)
        if gridM[1]:
            ax3.grid(axis='y',alpha=0.5, c='k', which=whichTicks)

#Helper function for resinv_plot
def __plot_pretty(inv_dict, x,z,v,im,cbarTicks,fig,ax, colMap='jet',cMin=None,cMax=None, gridFt=[False,False], gridM=[False,False], t='', primaryUnit='m', tight_layout=True, cBarOrient='vertical', cBarFormat ='%3.0f',cBarLabel ='Resistivity (ohm-m)', showPoints=False, norm=0, whichTicks='major', reverse=False):
    """Helper function for resinv_plot, parameters derived from there."""

    topoDF = inv_dict['topoDF']
    if cMin is None:
        cMin = inv_dict['resistModelDF']['Data'].min()
    if cMax is None:
        cMax = inv_dict['resistModelDF']['Data'].max()
    
    plt.rcParams["figure.dpi"] = 300  
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = plt.rcParams["xtick.top"] = True
    plt.sca(ax)
    

    vmax90 = np.percentile(v, 90)
    vmin2 = np.percentile(v, 2)
    vmax = v.max()
    vmin = v.min()
    minx = topoDF['xDist'].min()
    maxx = topoDF['xDist'].max()
    minz = min(z)
    maxz = max(z)

    #xlocsM = ax.get_xticks()
    if maxx>800:
        xlocsM = np.arange(minx,maxx+1,100)
        if max(xlocsM) < maxx:
            xlocsM = np.arange(minx,maxx+101,100)
    else:
        xlocsM = np.arange(minx,maxx+1,50)
        if max(xlocsM) < maxx:
            xlocsM = np.arange(minx,maxx+51,50)

    xlabelsM = [str(int(x)) for x in xlocsM]
    xlocsFt = np.uint16(xlocsM*3.2808399)
    xlabelsFt = [str(int(x)) for x in xlocsFt]

    minFtxLoc = xlocsFt[0]
    maxFtxLoc = xlocsFt[-1]
    xLabelsFtEven = np.arange(np.round(minFtxLoc,-1), np.round(maxFtxLoc,-1), 100)

    if np.round(maxFtxLoc,-1) < maxFtxLoc:
        xLabelsFtEven=np.insert(xLabelsFtEven, len(xLabelsFtEven),np.round(maxFtxLoc,-2))
        xLabelsFtEven=np.insert(xLabelsFtEven, len(xLabelsFtEven),np.round(maxFtxLoc,-2)+100)
    xLocs_FTinM = xLabelsFtEven / 3.2808399

    if np.ceil(maxz)>np.round(maxz,-1):
        zEnd = np.round(maxz,-1)+11
    else:
        zEnd = np.round(maxz,-1)+1
    
    ylocsM = np.arange(np.round(minz,-1),zEnd,10)
    ylabelsM = [str(x) for x in ylocsM]
    yLabelsFt = np.uint16(ylocsM*3.2808399)

    minFtyLabel = yLabelsFt.min()
    maxFtyLabel = yLabelsFt.max()
    yLabelsFtEven = np.arange(0, np.round(maxFtyLabel,-1), 20)
    

    if np.round(maxFtyLabel,-1) < maxFtyLabel:
        yLabelsFtEven=np.insert(yLabelsFtEven, len(yLabelsFtEven),yLabelsFtEven[-1]+20)
    yLocs_FTinM = yLabelsFtEven / 3.2808399

    yLimsM = [ylocsM[1],maxz+3]
    yLimsM = [minz,maxz+3]
    yLimsFt = [np.round(yLimsM[0]*3.2808399, 0), np.round(yLimsM[1]*3.2808399, 0)]

    ax.fill_between(topoDF['xDist'],topoDF['Elevation'],topoDF['Elevation']+10,color='w')
    ax.plot(topoDF['xDist'],topoDF['Elevation'],color='k',linewidth=1)
    ax.scatter(topoDF['xDist'],topoDF['Elevation'],marker='v',edgecolors='w',color='k',s=30)
    if showPoints:
        plt.scatter(x,z, c=v, marker='.', cmap='nipy_spectral', norm=norm)

    xLimsM = [minx,maxx]
    xLimsFt = [np.round(xLimsM[0]*3.2808399, 0), np.round(xLimsM[1]*3.2808399, 0)]
    if reverse:
        xLimsM.reverse()
        xLimsFt.reverse()
        #xlocsM=np.flip(xlocsM)
        #xLocs_FTinM=np.flip(xLocs_FTinM)


    if primaryUnit in ['m', 'meters', 'meter', 'metres', 'metre', 'metric']:
        pUnit = 'm'
        pUnitXLocs = xlocsM
        pUnitYLocs = ylocsM
        pUnitXLabels = xlabelsM
        pUnitYLabels = ylabelsM

        sUnit = 'ft'
        sUnitXLocs = xLocs_FTinM
        sUnitYLocs = yLocs_FTinM
        sUnitXLabels = xLabelsFtEven
        sUnitYLabels = yLabelsFtEven
        __label_plot(fig, ax, gridM, gridFt, whichTicks, primaryUnit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)
    
    elif primaryUnit in ['f', 'ft', 'feet', 'foot', 'US']:
        pUnit = 'ft'
        pUnitXLocs = xLocs_FTinM
        pUnitYLocs = yLocs_FTinM
        pUnitXLabels = xLabelsFtEven
        pUnitYLabels = yLabelsFtEven
        
        sUnit = 'm'
        sUnitXLocs = xlocsM
        sUnitYLocs = ylocsM
        sUnitXLabels = xlabelsM
        sUnitYLabels = ylabelsM
        __label_plot(fig, ax, gridM, gridFt, whichTicks, primaryUnit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)  
    
    else:
        pUnit = 'm'
        pUnitXLocs = xlocsM
        pUnitYLocs = ylocsM
        pUnitXLabels = xlabelsM
        pUnitYLabels = ylabelsM

        sUnit = 'ft'
        sUnitXLocs = xLocs_FTinM
        sUnitYLocs = yLocs_FTinM
        sUnitXLabels = xLabelsFtEven
        sUnitYLabels = yLabelsFtEven
        __label_plot(fig, ax, gridM, gridFt, whichTicks, primaryUnit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)
    
    if tight_layout:
        fig.tight_layout()
    else:
        plt.subplots_adjust(top=0.99, bottom=0.01, left=0.01, right=0.99, 
                    wspace=0, hspace=0)

    if cBarOrient=='horizontal':
        aspect=50
    else:
        aspect=25
    cbar = fig.colorbar(im, ax=ax,orientation=cBarOrient, aspect=aspect, extend='both',ticks=cbarTicks,format=cBarFormat)
    cbar.ax.tick_params(labelsize=12)#set_ticks(cbarTicks)
    cbar.set_label(label=cBarLabel, size=16)

    ax_h, ax_w = ax.bbox.height, ax.bbox.width
    axRatio = ax_w/ax_h
    aspRatio = (xLimsM[1]-xLimsM[0])/(yLimsM[1]-yLimsM[0])
    vertExag = abs(round(aspRatio/axRatio, 1))
    iterNo = inv_dict['iterationNo']
    iterInd = inv_dict['iterationInd']
    iterErr = inv_dict['iterationDF'].loc[iterInd, '%AbsError']
    ax.annotate('Iteration {}\nRMS Error: {}%'.format(iterNo, iterErr),xy=(0.02,0.02),xycoords='subfigure fraction', ha='left', va='bottom')
    ax.annotate('Vert.Exag: '+str(vertExag)+'x',xy=(0.98,0.02),xycoords='subfigure fraction', ha='right', va='bottom')
    fig.set_facecolor("w")

    return fig, ax
