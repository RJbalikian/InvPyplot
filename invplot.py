import csv
import pathlib 
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate


def autoplot(inv_file, iteration, return_dict=False, instrument='LS', auto_shift=True, **kwargs):
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
        If return_dict set to True, dictionary containing all input parameters and data generated along the way, including the output figures and axes

    """
    if isinstance(inv_file, pathlib.PurePath):
        pass
    else:
        inv_file = pathlib.Path(inv_file)

    inv_dict = ingest_inv(inv_file, verbose=False, instrument=instrument, show_iterations=False)
    inv_dict = read_inv_data(inv_file=inv_file, instrument=instrument, inv_dict=inv_dict)

    allIterList = [':', 'all']
    if type(iteration) is list:
        pass
    elif type(iteration) is int:
        iteration = [iteration]
    elif iteration.lower() in allIterList:
        iteration = inv_dict['iterationDF'].Iteration.tolist()

    resinv_params_list = ['inv_dict', 'colormap', 'cbar_format', 'cbar_label', 'cbar_orientation', 'cmin', 'cmax', 'auto_shift', 
                          'grid_ft', 'grid_m', 'title', 'norm_type', 'primary_unit', 'show_points','which_ticks', 
                          'figsize', 'dpi', 'reverse', 'tight_layout', 'save_fig', 'saveformat']
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
        inv_dict = get_resistivitiy_model(inv_file=inv_file, inv_dict=inv_dict, instrument=instrument)
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
def resinv_plot(inv_dict, title=None, primary_unit='m', colormap='nipy_spectral', grid_ft=[False,False], grid_m=[False,False], reverse=False, 
                cbar_orientation='horizontal', cmin=None, cmax=None, norm_type='log', cbar_format ='%3.0f', cbar_label='Resistivity (ohm-m)',
                auto_shift=True, show_points=False, which_ticks='major', figsize=None, dpi=None, tight_layout=True, save_fig=False, saveformat='png', imshow_kwargs={}, **kwargs):
    """Function to pull everything together and plot it nicely.

    It is recommended to use the autoplot function rather than resinv_plot directly, since using autoplot() incorporates all the setup needed to create the input dictionary keys/values correctly.

    Parameters
    ----------
    inv_dict : dict 
        Dictionary of inversion results generated from previous steps
    colormap : str, optional
        Colormap, any acceptable from matplotlib, by default 'nipy_spectral'
    cbar_format : str, optional
        Format string for colorbar tick labels, by default '%3.0f'
    cbar_label : str, optional
        Colorbar label, by default 'Resistivity (ohm-m)'
    cbar_orientation : str {'horizonta', 'vertical'}, optional
        Orientation of the colorbar, by default 'horizontal'
    cmin : float, optional
        Minimum of colorbar/colormap, by default None, which uses the minimum value of the dataset.
    cmax : float, optional
        Maximum of colorbar/colormap, by default None, which uses the maximum value of the dataset.
    grid_ft : list, optional
        Whether to show gridlines on the feet ticks, first position is x, second position is y, by default [False,False]
    grid_m : list, optional
        Whether to show gridlines on the meter tickes, first position is x, second posistion is y, by default [False,False]
    title : str, optional
        String to show as the title, if desired to set manually, by default None, which shows the filename as the title
    norm_type : str {'log', 'linear'}, optional
        Normalization type, by default 'log'. Determines whether matplotlib.colors.LogNorm or matplotlib.colors.Normalize is used for colormap.
    primary_unit : str {'m', 'ft'}, optional
        Whether to display meters or feet as primary unit (this determines which unit is larger on the axis and is on the left and top), by default 'm'
    show_points : bool, optional
        Whether to show the datapoints used for interpolation, by default False
    which_ticks : str {'major', 'minor', 'both'}, optional
        If grid_ft or grid_m has any True, this determines whether major, minor, or both gridlines are used; by default 'major'.
    figsize : tuple, optional
        Tuple (width, height) of the figsize, read into plt.rcParams['figure.figsize'], by default None.
    dpi : int or float, optional
        Resolution (dots per square inch) of final figure, read into plt.rcParams['figure.dpi'], by default None.
    reverse : bool, optional
        Whether to display the data in reverse (flipped along x) of what is read into from .inv file, by default False
    tight_layout : bool, optional
        If true, calls fig.tight_layout(). Otherwise, tries to maximize space on the figure using plt.subplots_adjust, by default True
    save_fig : bool, optional
        If False, will not save figure. Otherwise, calls plt.savefig() and the value of this parameter will be used as the output filepath, by default False.
    saveformat : str, optional
        Read into plt.savefig(format) paramater, by default 'png'.
    auto_shift : bool, default=True
        Whether to automatically shift the xDistances so the first one is at 0
    **imshow_kwargs
        Any extra specified keyword arguments are passed directly to matplotlib.pyplot.imshow

    Returns
    -------
    dict
        Returns existing inv_dict input, but with added keys of ['fig'] and ['ax'] containing a list of the fig and ax objects (list, since multiple iterations can be done at once)
    """
    #First get the data we'll need for the plot
    x = inv_dict['resistModelDF']['x'].copy().astype(np.float32)
    z = inv_dict['resistModelDF']['zElev'].copy().astype(np.float32)
    v = inv_dict['resistModelDF']['Data'].copy().astype(np.float32)

    if auto_shift:
        minXVal = min(x)
        x = x - minXVal
        inv_dict['topoDF']['xDist'] = inv_dict['topoDF']['xDist'] - minXVal

    #Setup plot parameters next
    if title is None:
        title = inv_dict['inv_file_Path'].stem
    
    if 'figure.dpi' not in list(inv_dict.keys()):
        inv_dict['figure.dpi'] = 250
    if 'figure.figsize' not in list(inv_dict.keys()):
        inv_dict['figure.figsize'] = (12,5)

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
    
    maxXDist = max(np.float64(inv_dict['electrodes']))

    if cmin is None:
        cmin = inv_dict['resistModelDF']['Data'].min()
    if cmax is None:
        cmax = inv_dict['resistModelDF']['Data'].max()
    
    #Convert data values to absolute values of float (need to remove abs()?)
    for i, val in enumerate(v):
        v[i] = abs(float(val))

    xi, zi = np.linspace(min(x), max(x), int(max(x))), np.linspace(min(z), max(z), int(max(abs(z))))
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

    vmax = cmax
    vmin = cmin

    #if cmax >= resistModelDF['Data'].max():
    #  vmax = vmax98
    #else:
    #  vmax = cmax
    #if cmin <= resistModelDF['Data'].min():
    #  vmin = vmin2
    #else:
    #  vmin = cmin
    #cbarTicks = np.arange(np.round(vmin,-1),np.round(vmax-1)+1,10) 

    arStep = np.round((vmax-vmin)/10,-1)
    cbarTicks = np.arange(np.round(vmin, -1), np.ceil(vmax/10)*10,arStep)

    #Get default values or kwargs, depending on if kwargs have been used
    if 'norm' in imshow_kwargs.keys():
        norm = imshow_kwargs['norm']
    else:
        if norm_type=='log':
            if vmin <= 0:
                vmin = 0.1
            norm = matplotlib.colors.LogNorm(vmin = vmin, vmax = vmax)
            #cbar_format = '%.1e'
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
        cmap=colormap   
        
    if 'interpolation' in imshow_kwargs.keys():
        interp = imshow_kwargs['interpolation']
        imshow_kwargs.pop('interpolation', None)
    else:
        interp='spline36'   

    imshow_kwargs.pop('imshow_kwargs', None)

    im = axes.imshow(vi, origin='lower',
                extent=extent,
                aspect=aspect,
                cmap =cmap,
                norm = norm,
                interpolation=interp, **imshow_kwargs)
    f, a = __plot_pretty(inv_dict, x,z,v,fig=fig,im=im,ax=axes,colormap=colormap,cmin=cmin,cmax=cmax, 
                       gridM=grid_m, gridFt=grid_ft, primary_unit=primary_unit, t=title, tight_layout=tight_layout,
                       cbarTicks=cbarTicks,cbar_format=cbar_format,cbar_label=cbar_label,cBarOrient=cbar_orientation,
                       show_points=show_points, norm=norm, which_ticks=which_ticks, reverse=reverse)

    plt.show()
    if save_fig is not False:
        f.savefig(save_fig, format=saveformat, facecolor='white')

    plt.close(f)
    return f, a

#Function to ingest inv file and find key information for later
def ingest_inv(inv_file, instrument='LS', verbose=True, show_iterations=True):
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

    global lsList
    global sasList

    lsList = ['terrameter ls', 'ls', 'terrameter', 'abem'] #Default, since it is current standard
    sasList= ['sas', '4000', 'terrameter sas4000', 'terrameter sas 4000']

    with open(str(inv_file)) as datafile: 
        #These may not be present without topography, initialize as None
        topoDataRow=None
        shiftMatrixRow=None
        electrodeCoordsRow=None
        noTopoPts=None

        #If not sensitivity or uncertainty calculated
        sensAndUncertainRow=None

        filereader = csv.reader(datafile)
        for row in enumerate(filereader):
            startLayer = 0
            endLayer = 0
            lay = -1

            global fileHeaderRows
            if instrument.lower() in sasList:
                fileHeaderRows = 5
                keyList=['Name', 'NomElectrodeSpacing', 'ArrayCode', 'NoDataPoints','DistanceType', 'FinalFlag']
            else:
                fileHeaderRows = 8
                keyList=['Name', 'NomElectrodeSpacing', 'ArrayCode', 'ProtocolCode', 'MeasureHeader', 'MeasureType', 'NoDataPoints','DistanceType','FinalFlag']

            if row[0] <= fileHeaderRows:
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
                #NOT CURRENTLY USED
                refResist = float(row[1][0].strip())
                continue
            
            if 'TOPOGRAPHICAL DATA' in str(row[1]):
                topoDataRow = row[0]
                continue
            if topoDataRow is not None:
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
    
    global use_topo
    if topoDataRow is None:
        use_topo=False
    else:
        use_topo=True
    return inv_dict

#Input Data
def read_inv_data(inv_file, inv_dict, instrument='LS'):
    """Function to do initial read of .inv file. 
    
    This data does not change with iteration, as in later function. This function should be readafter ingest_inv, using the output from that as inv_dict.
     
    

    Parameters
    ----------
    inv_file : str or pathlib.Path object
        Filepath of .inv file
    inv_dict : dict
        Dictionary (output from ingest_inv) containing information about where data is located in the .inv file
    instrument : str
        Which instrument was this data acquired with (changes a few of the read parameters)

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

    if instrument.lower() in sasList:
        startRow = 6
        inDF_cols=['xDist', 'aSpacing', 'nValue', 'Data']
        reSplitStr = ',\s+'
    elif instrument.lower() in lsList:
        startRow = 9
        inDF_cols=['NoElectrodes', 'A(x)', 'A(z)', 'B(x)', 'B(z)', 'M(x)', 'M(z)', 'N(x)', 'N(z)', 'Data']
        reSplitStr = '\s+'
    else:
        startRow = 9
        
        
    import csv
    with open(inv_file) as datafile: 
        filereader = csv.reader(datafile)
        start = 0
        inDataList = []
        for row in enumerate(filereader):
            if row[0] < startRow:
                continue
            elif row[0] < startRow+noDataPoints:
                if len(row[1]) == 1:
                    inDataList.append(re.sub(reSplitStr,' ',row[1][0]).split(' '))
                else:
                    for i, val in enumerate(row[1]):
                        row[1][i] = float(val.strip())
                    inDataList.append(row[1])
            else:
                break
    inDF = pd.DataFrame(inDataList)
    if startRow == 9:
        inDF.drop([0],inplace=True,axis=1)
        
    if inDF.shape[1] > 10:
        inDF = inDF.loc[:,0:10]
    inDF.replace(r'^\s*$', np.nan, regex=True, inplace=True)
    inDF.astype(np.float64)
    
    #Same older inv files have 3 columns only
    if inDF.shape[1] == 3:
        inDF.columns = ['xDist', 'aTimesN', 'Data']
        minASpacing = inDF['aTimesN'].min()   
        #minASpacing = inv_dict['fileHeaderDict']['NomElectrodeSpacing']
        
        inDF['aSpacing'] = minASpacing #Start with assumption that all are minASpacing
        
        for i in range(2,50):
            inDF['nValue'] = inDF['aTimesN'] / inDF['aSpacing']
            inDF['nTooBig'] = inDF['nValue'] >= 14 #Arbitrarily using value as the largest n value
            if inDF['nTooBig'].sum() == 0:
                break
            inDF.loc[inDF['nTooBig'] == True, 'aSpacing'] *= i #Wherever nTooBig, make aSpacing larger x2   
        inDF_cols = ['xDist', 'aTimesN', 'Data', 'aSpacing', 'nValue', 'nTooBig']
    inDF.columns = inDF_cols

    if instrument.lower() in sasList:   
        inDF['NoElectrodes'] = 4
        inDF['A(x)'] = inDF['xDist'] - inDF['aSpacing']/2
        inDF['A(z)'] = 0
        inDF['B(x)'] = inDF['A(x)'] + inDF['aSpacing']
        inDF['B(z)'] = 0        
        inDF['M(x)'] = inDF['A(x)'] + inDF['aSpacing'] + (inDF['aSpacing'] * inDF['nValue']) + inDF['aSpacing']
        inDF['M(z)'] = 0        
        inDF['N(x)'] = inDF['A(x)'] + inDF['aSpacing'] + (inDF['aSpacing'] * inDF['nValue'])
        inDF['N(z)'] = 0        

        inDF['xDist'] = pd.concat([inDF['xDist'], 
                                  inDF['xDist'] + inDF['aSpacing'],
                                  inDF['xDist'] + inDF['aSpacing'] + (inDF['aSpacing'] * inDF['nValue']) + inDF['aSpacing'],
                                  inDF['xDist'] + inDF['aSpacing'] + (inDF['aSpacing'] * inDF['nValue'])], ignore_index=True)
        inv_dict['xDists'] = pd.DataFrame(np.sort(inDF['xDist'].unique()), columns=['xDists'])

        inDF = inDF[['NoElectrodes', 'A(x)', 'A(z)', 'B(x)', 'B(z)', 'M(x)', 'M(z)', 'N(x)', 'N(z)', 'Data']]
    
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

    #shiftMatrixRow = None
    if use_topo: #inv_dict['shiftMatrixRow'] is not None:
        shiftMatrixRow = inv_dict['shiftMatrixRow']
        electrodeCoordsRow = inv_dict['electrodeCoordsRow']
        topoDataRow = inv_dict['topoDataRow']
        noTopoPts = inv_dict['noTopoPts']
    
    sensAndUncertainRow = inv_dict['sensAndUncertainRow']

    #Get Electrodes
    electrodes= pd.concat([invDF['A(x)'],invDF['B(x)'], invDF['M(x)'], invDF['B(x)'], invDF['N(x)']],ignore_index=True)
    electrodes.reset_index(inplace=True, drop=True)
    electrodes = electrodes.unique().astype(np.float32)
    inv_dict['electrodes'] = electrodes
    #inv_dict['xDists']  = electrodes+np.diff(electrodes)/2
    
    #ElectrodeCoordinates
    if use_topo:
        if shiftMatrixRow is not None:
            noModelElects = shiftMatrixRow-electrodeCoordsRow-1
        else:
            noModelElects = sensAndUncertainRow-electrodeCoordsRow-1
            
        electrodeCoordsDF = pd.read_table(inv_file,skiprows=electrodeCoordsRow,nrows=noModelElects, sep='\s+')
        electrodeCoordsDF.dropna(axis=1,inplace=True)
        electrodeCoordsDF.columns=['xDist','RelElevation']
        electrodeCoordsDF['ElectrodeNo'] = electrodeCoordsDF.index+1
    else:
        electrodeCoordsDF = pd.DataFrame({'xDist':electrodes})
        electrodeCoordsDF['RelElevation'] = 0.0
        electrodeCoordsDF['ElectrodeNo'] = electrodeCoordsDF.index+1       
    inv_dict['electrodeCoordsDF'] = electrodeCoordsDF

    #Topographical Data
    if use_topo:
        topoDF = pd.read_table(inv_file,skiprows=topoDataRow+2,nrows=noTopoPts, sep='\s+')
        topoDF.reset_index(inplace=True)
        topoDF.columns=['xDist','Elevation']
        topoDF['ElectrodeNo'] = topoDF.index+1
    else:
        topoDF = pd.DataFrame({'xDist':electrodes})
        topoDF['Elevation'] = 0
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
    if use_topo:
        if shiftMatrixRow is not None:
            shiftMatrixDF = pd.read_table(inv_file,skiprows=shiftMatrixRow+1,nrows=noModelElects, sep='\s+',header=None,index_col=0)
            shiftMatrixDF.dropna(axis=1,inplace=True)
            for c in shiftMatrixDF:
                shiftMatrixDF.rename(columns={c:'Layer'+str(int(c)-1)},inplace=True)
            inv_dict['shiftMatrixDF'] = shiftMatrixDF #Not currently using this
        else:
            #since shiftMatrixDF is not currently being used, maybe don't need to worry about this?
            pass

    #Sensitivity
    if sensAndUncertainRow is not None:
        sensDF = pd.read_table(inv_file,skiprows=sensAndUncertainRow+3,nrows=noModelBlocks, sep='\s+',header=None,index_col=0)
        sensDF.dropna(axis=1,inplace=True)
        sensDF.reset_index(inplace=True)
        sensDF.columns=['BlockNo','Sensitivity', '%ApproxUncertainty']
    else:
        sensDF = pd.DataFrame()
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
def get_resistivitiy_model(inv_file, inv_dict, instrument='LS'):
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
        
        currDF = pd.read_table(inv_file, skiprows=r[1]+1, nrows=noPtsInLyr, sep=',')
        #Clean up reading of table
        #Last line "x" value is blank, so make sure what is read in as x is actually data
        currDF.iloc[currDF.shape[0]-1,1] =  currDF.iloc[currDF.shape[0]-1,0]
        #Add in an "x" value for the last point (since we're making it tabular)
        currDF.iloc[currDF.shape[0]-1,0] = currDF.iloc[currDF.shape[0]-2,0]+1
        
        if instrument.lower() in sasList:
            currDF.columns = ['xDist','Data']
        else:
            currDF.columns=['ElectrodeNo','Data']
    
        currDF['zDepth'] = layerDepth

        for i in currDF.index:
            if instrument.lower() in lsList:
                lowerElecNo = currDF.loc[i,'ElectrodeNo']#-1
                elecInd = electrodeCoordsDF.loc[electrodeCoordsDF['ElectrodeNo']==lowerElecNo].index.values[0]
                currDF.loc[i,'x'] = (electrodeCoordsDF.loc[elecInd,'xDist'] + electrodeCoordsDF.loc[elecInd+1,'xDist'])/2
            else:
                currDF.loc[i,'x'] = currDF.loc[i,'xDist']
                
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
def __label_plot(fig, ax, gridM, gridFt, which_ticks, pUnit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims, yLims, t):
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
           ax2.grid(axis='x',alpha=0.5, c='k', which=which_ticks)
        if gridM[1]:
            ax.grid(axis='y',alpha=0.5, c='k', which=which_ticks)
        if gridFt[0]:
            ax.grid(axis='x',alpha=0.5, c='k', which=which_ticks)
        if gridFt[1]:
            ax3.grid(axis='y',alpha=0.5, c='k', which=which_ticks)
    else:
        if gridFt[0]:
            ax2.grid(axis='x',alpha=0.5, c='k', which=which_ticks)
        if gridFt[1]:
            ax.grid(axis='y',alpha=0.5, c='k', which=which_ticks)
        if gridM[0]:
            ax.grid(axis='x',alpha=0.5, c='k', which=which_ticks)
        if gridM[1]:
            ax3.grid(axis='y',alpha=0.5, c='k', which=which_ticks)

#Helper function for resinv_plot
def __plot_pretty(inv_dict, x,z,v,im,cbarTicks,fig,ax, colormap='nipy_spectral',cmin=None,cmax=None, gridFt=[False,False], gridM=[False,False], t='', primary_unit='m', tight_layout=True, cBarOrient='vertical', cbar_format ='%3.0f',cbar_label ='Resistivity (ohm-m)', show_points=False, norm=0, which_ticks='major', reverse=False):
    """Helper function for resinv_plot, parameters derived from there."""

    topoDF = inv_dict['topoDF']
    if cmin is None:
        cmin = inv_dict['resistModelDF']['Data'].min()
    if cmax is None:
        cmax = inv_dict['resistModelDF']['Data'].max()
    
    plt.rcParams["figure.dpi"] = 300  
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = plt.rcParams["xtick.top"] = True
    plt.sca(ax)
    

    vmax90 = np.percentile(v, 90)
    vmin2 = np.percentile(v, 2)
    vmax = v.max()
    vmin = v.min()
    minx = min(x)#topoDF['xDist'].min()
    maxx = max(x)#topoDF['xDist'].max()
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
    if show_points:
        plt.scatter(x,z, c=v, marker='.', cmap=colormap, norm=norm)

    xLimsM = [minx,maxx]
    xLimsFt = [np.round(xLimsM[0]*3.2808399, 0), np.round(xLimsM[1]*3.2808399, 0)]
    if reverse:
        xLimsM.reverse()
        xLimsFt.reverse()
        #xlocsM=np.flip(xlocsM)
        #xLocs_FTinM=np.flip(xLocs_FTinM)


    if primary_unit in ['m', 'meters', 'meter', 'metres', 'metre', 'metric']:
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
        __label_plot(fig, ax, gridM, gridFt, which_ticks, primary_unit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)
    
    elif primary_unit in ['f', 'ft', 'feet', 'foot', 'US']:
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
        __label_plot(fig, ax, gridM, gridFt, which_ticks, primary_unit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)  
    
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
        __label_plot(fig, ax, gridM, gridFt, which_ticks, primary_unit, pUnitXLocs, pUnitYLocs, pUnitXLabels,pUnitYLabels, sUnit, sUnitXLocs, sUnitYLocs, sUnitXLabels, sUnitYLabels, xLims=xLimsM, yLims=yLimsM, t=t)
    
    if tight_layout:
        fig.tight_layout()
    else:
        plt.subplots_adjust(top=0.99, bottom=0.01, left=0.01, right=0.99, 
                    wspace=0, hspace=0)

    if cBarOrient=='horizontal':
        aspect=50
    else:
        aspect=25
    cbar = fig.colorbar(im, ax=ax,orientation=cBarOrient, aspect=aspect, extend='both',ticks=cbarTicks,format=cbar_format)
    cbar.ax.tick_params(labelsize=12)#set_ticks(cbarTicks)
    cbar.set_label(label=cbar_label, size=16)

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
