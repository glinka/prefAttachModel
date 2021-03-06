from multiprocessing import Pool
import getGraphProps as gGP
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter, FormatStrFormatter
import numpy as np
import os
import util_fns as uf

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from pca import pca

import pygraphviz as pgv

# pgv.Draw

# from mpltools import style
# from mpltools import layout

# style.use('ggplot')

# global variables to the (hopefully temporary) rescue of aligning colornorms
# PLEASE DELETE ASAP

def gamma_bullshit():
    """gamma function investigatoin"""
    f = lambda x, k: np.power(x, k-1)*np.exp(-x)/scipy.special.gamma(k)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(0, 120, 100), f(np.linspace(0, 120, 100), 10))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(0, 120, 100), f(np.linspace(0, 120, 100), 20))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(0, 120, 100), f(np.linspace(0, 120, 100), 100))
    plt.show()
    
def thin_array(array, frac_to_keep=0.5, new_npts=None):
    #
    # !!! shape of array must be (x, y), cannot be (x,) !!!
    #
    # will keep 100*frac_to_keep% of the data 
    # in each column of array, evenly spaced
    # thus, array has form col1 col2 col3 .. colN
    # and returns col1 col2...colN reduced to frac_to_keep of
    # original length
    # no thinning is possible for frac_to_keep > 0.5,
    # at least in this arrangement
    if new_npts is None:
        if frac_to_keep > 0.5:
            return array
        else:
            npts = array.shape[0]
            new_npts = int(frac_to_keep*npts)
            spacing = npts/new_npts
            ncols = array.shape[1]
            thinned_array = np.zeros((new_npts, ncols))
    else:
        npts = array.shape[0]
        if new_npts > npts/2:
            return array
        else:
            spacing = npts/new_npts
            ncols = array.shape[1]
            thinned_array = np.zeros((new_npts, ncols))
    for i in range(new_npts):
        thinned_array[i,:] = array[spacing*i, :]
    return thinned_array

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    f = open(path_to_file, "r")
    params_str = f.readline()
    params = get_header_data(params_str)
    f.close()
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    print 'done importing data from', filename
    return data, params

def get_header_data(header_str):
    BEGIN = 0
    comma = 1
    #create dict from header, based on key=value format in csv
    params = {}
    while comma > 0:
        equals = header_str.find("=")
        comma = header_str.find(",")
        params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
        header_str = header_str[comma+1:]
    params[header_str[BEGIN:equals]] = float(header_str[equals+1:comma])
    #make integer, may not work especially well
    for key in params:
        if(params[key] % 1 == 0):
            params[key] = int(params[key])
    return params

def genData(fileName):
    pathToFile = os.path.realpath(fileName)
    f = open(pathToFile, "r")
    paramstr = f.readline()
    params = {}
    f.close()
    begin = 0
    comma = 1
    #create dict from header, based on key=value format in csv
    while comma > 0:
        comma = paramstr.find(",")
        equals = paramstr.find("=")
        params[paramstr[begin:equals]] = float(paramstr[equals+1:comma])
        paramstr = paramstr[comma+1:]
    params[paramstr[begin:equals]] = float(paramstr[equals+1:comma])
    for key in params:
        if(params[key] % 1 == 0):
            params[key] = int(params[key])
    data = np.genfromtxt(pathToFile, delimiter=",", skip_header=1)
    print params
    return params, data

def makeFolder(baseName):
    folderCount = 0
    newFolder = baseName + str(folderCount) + '/'
    while os.path.exists(newFolder):
        folderCount = folderCount + 1
        newFolder = baseName + str(folderCount) + '/'
    os.mkdir(newFolder)
    return newFolder

def genFileName(baseName, params, uniqueID=''):
    fileName = baseName
    for key in params.keys():
        fileName = fileName + '_' + key + '_' + str(params[key])
    return fileName + '_' + uniqueID

def animateContour(data, params, cmap='Paired', fps=10, bitrate=14400, containerType='.mkv'):
    import matplotlib.animation as animation
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure(facecolor='w')
    plt.text(n/2,1.05*params['n'],'Evolution of adjacency matrix in preferential attachment model', ha='center', fontsize=16)
    newFolder = makeFolder('contourPlots')
    for i in range(nData):
        plt.clf()
        plt.pcolormesh(data[(i)*n:(i+1)*n,:n], cmap=cmap)
        plt.draw()
        fileName = genFileName('contour', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

def animateRawData(data, params, fps=10, bitrate=14400, containerType='.mkv'):
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import Axes3D
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure(facecolor='w')
    plt.figtext(0.5, 0.92 ,'Evolution of adjacency matrix in preferential attachment model', ha='center', fontsize=16)
    spAxes = [fig.add_subplot(i, projection='3d') for i in range(221, 225)]
    spAxes[0].view_init(-2.0, 45.0)
    spAxes[1].view_init(-2, 135)
    spAxes[2].view_init(45, 225)
    #find max degree to set z limits:
    maxDeg = np.amax(data[:,:])
    for ax in spAxes:
        ax.set_xlim3d(left=0, right=n)
        ax.set_ylim3d(bottom=0, top=n)
        ax.set_zlim3d(bottom=0, top=maxDeg)
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    newFolder = makeFolder('animations')
    for i in range(nData):
        for ax in spAxes:
            ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet')
        plt.draw()
        fileName = genFileName('rawData3d', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

def animateReconstruction(data, params, fps=10, bitrate=1600, containerType='.mkv'):
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import Axes3D
    n = params['n']
    nData = data.shape[0]/n
    fig = plt.figure(facecolor='w')
    ax1 = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212, projection='3d')
    ax1.grid(b=False)
    ax2.grid(b=False)
    #find max degree to set z limits:
    maxDeg = np.amax(data[:,:])
    ax1.set_xlim3d(left=0, right=n)
    ax1.set_ylim3d(bottom=0, top=n)
    ax1.set_zlim3d(bottom=0, top=maxDeg)
    ax2.set_xlim3d(left=0, right=n)
    ax2.set_ylim3d(bottom=0, top=n)
    ax2.set_zlim3d(bottom=0, top=maxDeg)
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    newFolder = makeFolder('reconstruction')
    print newFolder
    fileName = ""
    for i in range(nData):
        z = data[(i)*n:(i+1)*n,:n]
        zRecon = gGP.getEigenReconstruction(z)
        ax1.scatter(xgrid, ygrid, z, c=z, cmap='jet')
        ax2.scatter(xgrid, ygrid, zRecon, c=zRecon, cmap='RdBu')
        fileName = genFileName('rawDataAndReconstruction', params, str(i))
        plt.draw()
        plt.savefig(newFolder + fileName + '.png')
        ax1.cla()
        ax2.cla()
        ax1.grid(b=False)
        ax2.grid(b=False)
        ax1.set_xlim3d(left=0, right=n)
        ax1.set_ylim3d(bottom=0, top=n)
        ax1.set_zlim3d(bottom=0, top=maxDeg)
        ax2.set_xlim3d(left=0, right=n)
        ax2.set_ylim3d(bottom=0, top=n)
        ax2.set_zlim3d(bottom=0, top=maxDeg)
        print 1.0*(i+1)/nData
    makeAnimation(fileName, newFolder)

def makeAnimation(inputFilename, inputFolder, fps=50, bitrate=3000000, containerType='.mkv'):
    from subprocess import call
    us1 = 1
    us2 = 0
    while us1 > 0:
        us2 = us1
        us1 = inputFilename.find("_", us2+1)
    fileNamebase = inputFilename[:us2+1]
    inputFiles =  fileNamebase + "%d.png"
    outputFilename = fileNamebase + containerType
    os.chdir(os.path.realpath(inputFolder))
    call(["ffmpeg", "-i", inputFiles, "-r", str(fps), "-b", str(bitrate), outputFilename])

def plotCRecon(data, params):
    from mpl_toolkits.mplot3d import Axes3D
    n = params['n']
    nData = data.shape[0]/(2*n)
    fig = plt.figure(facecolor='w')
    ax1 = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212, projection='3d')
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    maxDeg = np.amax(data)
    ax1.set_xlim(left=0, right=n)
    ax1.set_ylim(bottom=0, top=n)
    ax1.set_zlim(bottom=0, top=maxDeg)
    ax2.set_xlim(left=0, right=n)
    ax2.set_ylim(bottom=0, top=n)
    ax2.set_zlim(bottom=0, top=maxDeg)
    newFolder = makeFolder('CRecon')
    print '--> saving', nData, 'images in', newFolder
    fileName = ''
    for i in range(nData):
        preRecon = data[n*2*i:n*(2*i+1),:]
        postRecon = data[n*(2*i+1):n*(2*i+2),:]
        ax1.scatter(xgrid, ygrid, preRecon, c=preRecon, cmap='jet')
        ax2.scatter(xgrid, ygrid, postRecon, c=postRecon, cmap='jet')
        plt.draw()
        fileName = genFileName('CRecon', params, uniqueID=str(i))
        plt.savefig(newFolder + fileName + '.png')
        ax1.cla()
        ax2.cla()
        ax1.grid(b=False)
        ax2.grid(b=False)
        ax1.set_xlim3d(left=0, right=n)
        ax1.set_ylim3d(bottom=0, top=n)
        ax1.set_zlim3d(bottom=0, top=maxDeg)
        ax2.set_xlim3d(left=0, right=n)
        ax2.set_ylim3d(bottom=0, top=n)
        ax2.set_zlim3d(bottom=0, top=maxDeg)
    makeAnimation(fileName, newFolder)

def compare_recon(data_list, params):
    # expects a bunch of projData csv files in which
    # data is layered in before_proj, after_proj
    # adjacency matrices, the degree distributions of which
    # are subtracted and plotted to give a sense of the errors
    # in the reconstruction process
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    n = params['n']
    ndata = len(data_list)
    print ndata
    npts = (data_list[0]).shape[0]/(2*n)
    err = np.empty(npts)
    max_err = np.empty(npts)
    avg_err = np.empty(npts)
    # avg the data, then plot
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    colornorm = colors.Normalize(vmin=0, vmax=npts-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    for i in range(npts):
        preRecon = np.zeros((n,n))
        postRecon = np.zeros((n,n))
        for k in range(ndata):
            preRecon = preRecon + data_list[k][n*2*i:n*(2*i+1),:]
            postRecon = postRecon + data_list[k][n*(2*i+1):n*(2*i+2),:]
        preRecon = preRecon / float(ndata)
        postRecon = postRecon / float(ndata)
        pre_recon_degs = np.sum(preRecon, 0)
        post_recon_degs = np.sum(postRecon, 0)
        if i == 0 or i == npts-1:
            ax.plot(range(n), np.sort(pre_recon_degs) -  np.sort(post_recon_degs), c=colormap.to_rgba(float(i)), lw=5)
        else:
            ax.plot(range(n),  np.sort(pre_recon_degs) -  np.sort(post_recon_degs), c=colormap.to_rgba(float(i)))
        pre_sort = np.argsort(pre_recon_degs)
        post_sort = np.argsort(post_recon_degs)
        preRecon = preRecon[pre_sort, :]
        preRecon = preRecon[:, pre_sort]
        postRecon = postRecon[post_sort, :]
        postRecon = postRecon[:, post_sort]
        err[i] = np.linalg.norm(preRecon - postRecon)
        max_err[i] = np.max(preRecon - postRecon)
        avg_err[i] = np.average(preRecon - postRecon)
    ax.set_xlabel('index', fontsize=24)
    ax.set_ylabel('degree difference (actual - reconstruction)', fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=24)
    plt.show()
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.plot(range(npts), err, c='c', label='cumulative')
    ax.plot(range(npts), avg_err, c='g', label='average')
    ax.plot(range(npts), max_err, c='r', label='maximum single-degree')
    ax.set_ylim(bottom=0)
    ax.set_xlabel('approximate time')
    ax.set_ylabel('reconstruction error')
    print '--> saving images in recon_comp'
    plt.savefig('./recon_comp/recon_comp.png')
    

def plotEigVectRecon(data, params):
    n = params['n']
    nData = data.shape[0]/(2)
    fig = plt.figure(facecolor='w')
    ax1 = fig.add_subplot(111)
    yMin = np.amin(data)
    yMax = np.amax(data)
    #ax1.set_ylim((yMin, yMax))
    xData = np.linspace(1, n, n)
    folderName = makeFolder("eigVectRecon")
    print folderName
    for i in range(nData):
        ax1.set_xlabel('index')
        ax1.set_ylabel('vector value')
        # ax1.plot(xData, data[2*i*n: n*(2*i+1)])
        # ax1.plot(xData, data[n*(2*i+1): 2*n*(i+1)], c='g')
        ax1.hold(True)
        ax1.plot(xData, data[i,:], c='k')
        ax1.plot(xData, data[nData+i,:], c='g')
        plt.savefig(folderName + genFileName("eigVectRecon", params, uniqueID=str(i)) + '.png')
        ax1.cla()

def plotFittedData(data, params, fns):
    from mpl_toolkits.mplot3d import Axes3D
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure(facecolor='w')
    fig.hold(True)
    spAxes = [fig.add_subplot(i, projection='3d') for i in range(221, 225)]
    spAxes[0].view_init(-2.0, 45.0)
    spAxes[1].view_init(-2, 135)
    spAxes[2].view_init(45, 225)
    xyScaling = 100.0
    xgrid = xgrid/xyScaling
    ygrid = ygrid/xyScaling
    #find max degree for axis limits
    maxDeg = np.amax(data[:,:])
    for ax in spAxes:
        ax.set_xlim3d(left=0, right=n/xyScaling)
        ax.set_ylim3d(bottom=0, top=n/xyScaling)
        ax.set_zlim3d(bottom=0, top=maxDeg)
    #don't plot each f(x,y) in wireframe
    stride = 10
    newFolder = makeFolder('fittedPlots')
    for i in range(nData):
        Z = np.transpose(data[(i)*n:(i+1)*n,:n])
        fn = gGP.fitXYFunction(xgrid, ygrid, Z, fns)
        for ax in spAxes:
            ax.scatter(xgrid, ygrid, Z, c=Z, cmap='jet', alpha=0.1)
            ax.plot_wireframe(xgrid, ygrid, fn(xgrid, ygrid), rstride=stride, cstride=stride)
        plt.draw()
        fileName = genFileName('fittedPlot', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

def animateVector(data, params, fn, fps=10, bitrate=14400, containerType='.mkv'):
    n = params['n']
    nData = data.shape[0]/n
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.set_xlabel('Vector index')
    ax.set_ylabel('Vector value')
    newFolder = makeFolder('vector')
    print newFolder
    fileName = ""
    #find y upper and lower limits
    yMin = np.min([fn(data[(i)*n:(i+1)*n,:n]) for i in range(nData)])
    yMax = np.max([fn(data[(i)*n:(i+1)*n,:n]) for i in range(nData)])
    for i in range(nData):    
        ax.cla()
        ax.set_xlabel('Vector index')
        ax.set_ylabel('Vector value')
        ax.set_ylim((yMin, yMax))
        ax.plot(np.linspace(1,n,n), fn(data[(i)*n:(i+1)*n,:n]), marker='o', c=[1,0.5,0.5])
        fileName = genFileName('eigVals', params, str(i))
        plt.savefig(newFolder + fileName + '.png')
    makeAnimation(fileName, newFolder)

def comp_eigvect_recon(coeffs, params, plot_name=""):
    n = params['n']
    proj_step = params['proj_step']
    wait = params['off_manifold_wait']
    nmicrosteps = params['nms']
    collect_interval = params['collection_interval']
    nsaves_per_proj = (nmicrosteps - wait)/collect_interval + 2
    nprojs = (data.shape[0] - 1)/nsaves_per_proj
    print nprojs, nsaves_per_proj
    nvects = data.shape[1]
    nyvects = nvects - 1
    if plot_name is not "":
        plot_name = plot_name + "_"
    line = np.arange(n)
    for i in range(nprojs):
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        for j in range(nsaves_per_proj):
            recon = np.zeros(n)
            for k in range(nyvects-1, -1, -1):
                recon = recon*line + coeffs[i*nsaves_per_proj + j][k]
            if j == nsaves_per_proj-1:
                ax.plot(range(n), recon, c='r')
            else:
                ax.plot(range(n), recon, c='b')
        ax.set_title('reconstructions over time')
        plt.savefig("eigvect_recon/" + plot_name + "eigvect_recon" + str(i) + ".png")

def animate_eigvals(eigvals, params, fps=10, bitrate=14400, containerType='.mkv'):
    n = params['n']
    nvects = eigvals.shape[0]
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.set_xlabel('Index')
    ax.set_ylabel('Eigenvalue')
    newFolder = makeFolder('eigenvalues')
    print 'saved in', newFolder
    fileName = ""
    #find y upper and lower limits
    yMin = np.min(eigvals)
    yMax = np.max(eigvals)
    for i in range(nvects):    
        ax.cla()
        ax.set_xlabel('Index')
        ax.set_ylabel('Eigenvalue')
        # ax.set_ylim((yMin, yMax))
        ax.plot(np.arange(n), np.sort(np.log(np.abs(data[i, :]))), marker='o', c=[1,0.5,0.5])
        fileName = genFileName('eigvals', params, str(i))
        plt.savefig(newFolder + fileName + '.png')
    makeAnimation(fileName, newFolder)

def scalarEvolution(data, params, fn):
    nData = params['nSteps']/params['dataInterval']
    stepSize = params['dataInterval']
    n = params['n']
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    newFolder = makeFolder('vector')
    fileName = ""
    yData = np.array([fn(data[(i)*n:(i+1)*n,:n]) for i in range(nData)])
    xData = np.array([i*stepSize for i in range(nData)])
    ax.plot(xData, yData, color='g', marker='o')
    fileName = genFileName('scalarEvo', params, str(i))
    newFolder = makeFolder('scalarEvo')
    plt.savefig(newFolder + fileName + ".png")
    

def plotReconstruction((x, y, z, zlim, xylim, folder, fileName)):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(facecolor='w')
    ax1 = fig.add_subplot(211, projection='3d')
    ax2 = fig.add_subplot(212, projection='3d')
    ax1.grid(b=False)
    ax2.grid(b=False)
    ax1.set_xlim3d(left=0, right=xylim)
    ax1.set_ylim3d(bottom=0, top=xylim)
    ax1.set_zlim3d(bottom=0, top=zlim)
    ax2.set_xlim3d(left=0, right=xylim)
    ax2.set_ylim3d(bottom=0, top=xylim)
    ax2.set_zlim3d(bottom=0, top=zlim)
    ax1.scatter(x, y, z, c=z, cmap='jet')
    ax2.scatter(x, y, gGP.getEigenReconstruction(z), c=z, cmap = 'jet')
    plt.savefig(folder + fileName + ".png")
    plt.close(fig)

def compareProjection(fullData, fullParams, cpiData, cpiParams):
    #3d plot of full simulation and cpi, degree evolution vs time
    from mpl_toolkits.mplot3d import Axes3D
    n = fullParams['n']
    nFullData = fullData.shape[0]/n
    nCPIData = cpiData.shape[0]/n
    #full collection interval
    fullCI = fullParams['dataInterval']
    #cpi collection interval
    cpiCI = cpiParams['collectInterval']
    cpiNMS = cpiParams['nMicroSteps']
    cpiOMW = cpiParams['offManifoldWait']
    cpiProjStep = cpiParams['projStep']
    #on manifold steps, with any luck is an integer
    cpiOMS = (cpiNMS - cpiOMW)/cpiCI
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlabel('degree')
    ax.set_ylabel('steps')
    ax.set_xlabel('vertex index')
    xData = np.linspace(1,n,n)
    for i in range(nFullData):
        time = i*fullCI*np.ones(n)
        ax.plot(xData, time, gGP.getDegrees(fullData[i*n:(i+1)*n,:]), c='r', alpha=0.1)
    t = cpiOMW
    onManifoldCount = 0
    for i in range(nCPIData):
        time = t*np.ones(n)
        if onManifoldCount < cpiOMS:
            t = t + cpiCI
            onManifoldCount = onManifoldCount + 1
        else:
            t = t + cpiProjStep + cpiOMW
            onManifoldCount = 0
        ax.plot(xData, time, gGP.getDegrees(cpiData[i*n:(i+1)*n,:]), c='b', alpha=0.1)
    newFolder = makeFolder('compProj')
    plt.savefig(newFolder + 'compProj.png')

def makeDegSurface(data, params, fn):
    from mpl_toolkits.mplot3d import Axes3D
    n = params['n']
    nData = data.shape[0]/n
    fig = plt.figure(facecolor='w')
    fig2 = plt.figure(facecolor='w')
    ax2 = fig2.add_subplot(111, projection='3d')
    fig.hold(True)
    #spAxes = [fig.add_subplot(i, projection='3d') for i in range(211, 213)]
    spAxes = [fig.add_subplot(111, projection='3d')]
    spAxes[0].view_init(30, -135)
    ax2.view_init(30, -135)
    maxZ = np.max([np.max(fn(data[i*n:(i+1)*n,:])) for i in range(nData)])
    ci = params['dataInterval']
    tSpan = nData*ci
    nNCubedSteps = tSpan/np.power(n, 3)
    nNSqSteps = np.power(n, 2)/ci
    #set various axis properties
    #[ax.grid(b=False) for ax in spAxes]
    [ax.set_xlim((0, n)) for ax in spAxes]
    [ax.set_xticklabels([str(n)]) for ax in spAxes]
    [ax.set_xticks([n]) for ax in spAxes]
    [ax.set_xlabel('vertex index') for ax in spAxes]
    [ax.set_ylim((1, tSpan)) for ax in spAxes]
    [ax.set_yticks([str(i) for i in (np.power(n, 3)*np.arange(nNCubedSteps))]) for ax in spAxes]
    [ax.set_yticklabels([str(i) for i in (np.power(n, 3)*np.arange(nNCubedSteps))]) for ax in spAxes]
    [ax.set_ylabel('simulation step') for ax in spAxes]
    ax2.set_xlabel('vertex index')
    ax2.set_ylabel('simulation step')
    ax2.set_zlabel('degree')
    # [ax.set_zlim((0, maxZ)) for ax in spAxes]
    # [ax.set_zticklabels([str(maxZ)]) for ax in spAxes]
    # [ax.set_zticks([maxZ]) for ax in spAxes]
    [ax.set_zlabel('degree') for ax in spAxes]
    # [ax.set_zscale('log') for ax in spAxes]
    xData = np.linspace(1, n, n)
    for i in range(nData):
        ys = ci*(i+1)*np.ones(n)
        color = 'b'
        if (i+1)*params['dataInterval'] % (np.power(params['m'], 3)) == 0:
            color = 'r'
        if((i+1)*ci % np.power(n, 3) == 0):
            [ax.plot(xData, ys, np.log(1+fn(data[i*n:(i+1)*n])), c=color, alpha=0.5) for ax in spAxes]
        if(i*ci < np.power(n, 3)):
            ax2.plot(xData, ys, np.log(1+fn(data[i*n:(i+1)*n])), c=color, alpha=0.5)
    newFolder =makeFolder('vectorSurface')
    fig.savefig(newFolder + 'degreeSurfacen3.png')
    fig2.savefig(newFolder + 'degreeSurfacen2.png')

def plot_coeffs(times, coeffs_list, plot_name=""):
    coeffs = np.average(np.array(coeffs_list), 0)
    n = coeffs.shape[0]
    ncoeffs = coeffs.shape[1]
    if plot_name is not "":
        plot_name = plot_name + "_"
    for i in range(ncoeffs):
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        ax.scatter(times, coeffs[:,i], lw=0)
        maxcoeff = np.amax(coeffs[:,i])
        mincoeff = np.amin(coeffs[:,i])
        ax.set_ylim((mincoeff - 0.1*(maxcoeff - mincoeff), maxcoeff + 0.1*(maxcoeff - mincoeff)))
        plt.savefig("coeffs/" + plot_name + "coeff" + str(i) + ".png")

def plot_coeffs_fitting(times, coeffs_list, fit, plot_name=""):
    coeffs = np.average(np.array(coeffs_list), 0)
    n = coeffs.shape[0]
    ncoeffs = coeffs.shape[1]
    nfitcoeffs = fit.shape[1]
    if plot_name is not "":
        plot_name = plot_name + "_"
    for i in range(ncoeffs):
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        ax.scatter(times, coeffs[:,i], lw=0, c='b')
        evals = np.zeros(n)
        for j in range(nfitcoeffs-1, 0, -1):
            evals = 10.0*((np.linspace(10000, 200000, n) - 10000)/80000)*(evals + fit[i,j])
        evals = evals + np.ones(n)*fit[i,0]
        ax.plot(np.linspace(10000, 200000, n), evals, color='g')
        maxcoeff = np.amax(coeffs[:,i])
        mincoeff = np.amin(coeffs[:,i])
        ax.set_ylim((mincoeff - 0.1*(maxcoeff - mincoeff), maxcoeff + 0.1*(maxcoeff - mincoeff)))
        ax.set_xlim((0,2*200000))
        plt.savefig("coeffs/" + plot_name + "coeff" + str(i) + ".png")

def plot_vectors_tc(data, params, plot_name=""):
    """ Assumes data is arranged as:

    vector1   vector2   vector3 ... vectorN     times
    data[0]   data[1]   data[2]     data[N-2]   data[N-1]

    and plots vectors against times
    """
    proj_step = params['proj_step']
    wait = params['off_manifold_wait']
    nmicrosteps = params['nms']
    collect_interval = params['collection_interval']
    nsaves_per_proj = (nmicrosteps - wait)/collect_interval + 2
    nprojs = (data.shape[0] - 1)/nsaves_per_proj
    print nprojs, nsaves_per_proj
    nvects = data.shape[1]
    nyvects = nvects - 1
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    if plot_name is not "":
        plot_name = plot_name + "_"
    for i in range(nyvects):
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        for j in range(nprojs):
            ax.scatter(data[j*nsaves_per_proj:(j+1)*nsaves_per_proj,nyvects], data[j*nsaves_per_proj:(j+1)*nsaves_per_proj,i], c=(np.sin(float(i)/nyvects), np.cos(float(1-i)/nyvects), 1-float(i)/nyvects), label="coeff: " + str(i+1), lw=0)
        ax.scatter(np.arange(1, nprojs+1)*(nmicrosteps + proj_step), data[np.arange(nprojs)*nsaves_per_proj + nsaves_per_proj - 1, i], lw=0, s=30, c='r')
        ax.set_xlabel('simulation step')
        ax.set_ylabel('coefficient value')
        ax.set_xlim(left=0)
        plt.savefig("coeffs/" + plot_name + "coeff" + str(i) + ".png")
    #ax.legend(loc=6)

def plot_degree_surface(degs, times, sort=True, title='', zlabel=None, ax=None, FONTSIZE=48, zlim=None, colornorm=None):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar

    LABELSIZE = 0.65*FONTSIZE
    LEGENDSIZE = 0.5*FONTSIZE
    ALPHA = 0.8

    new_npts = 50
    times.shape = (times.shape[0], 1)
    times = thin_array(times, new_npts=new_npts)
    degs = thin_array(degs, new_npts=new_npts)

    if sort:
        sorted_degs = np.sort(degs)
    else:
        sorted_degs = degs
    n = degs.shape[1]
    indices = np.linspace(1, n, n)
    maxdeg = np.amax(sorted_degs)
    mindeg = np.amin(sorted_degs)
    # color by deg (z val), not vertex (x val)
    # colornorm = colors.Normalize(vmin=0, vmax=n-1)
    if colornorm is None:
        colornorm = colors.Normalize(vmin=mindeg, vmax=maxdeg)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    ntimes = times.shape[0]

    show = False
    if ax is None:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111, projection='3d')
        show = True

    npts = times.shape[0]
    for v in range(n):
        ax.scatter(v*np.ones(npts), times, sorted_degs[:,v], linewidths=0, c='b')#colormap.to_rgba(1.0*sorted_degs[:,v])) # *v))
    ax.set_xlim((0,n))
    ax.set_xticklabels([str(int(i)) for i in np.linspace(1, n, 6)])
    ax.set_ylim(bottom=0)
    if zlim is None:
        zlim = (mindeg-(maxdeg-mindeg)*0.1, maxdeg+(maxdeg-mindeg)*0.1)
        ax.set_zlim(zlim)
    else:
        ax.set_zlim(zlim)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.set_xlabel('\nvertex', fontsize=FONTSIZE)
    ax.set_ylabel('\nstep', fontsize=FONTSIZE)
    if zlabel is None:
        zlabel = 'degree'
    ax.set_zlabel(zlabel, fontsize=FONTSIZE)
    ax.set_title(title, fontsize=FONTSIZE)
    # increase size of axes scaling number and add padding
    # scale_pos = ax.yaxis.get_children()[1].get_position()
    # print ax.yaxis.get_children()[1]
    # ax.yaxis.get_children()[1].set_text("\n" + str(scale_txt))
    ax.yaxis.get_children()[1].set_size(0)
    # SUPER DAMN STUPID, BUT I DON'T KNOW HOW TO PROGRAMMATICALLY FIND THE EXPONENTIAL SCALE VALUE
    # DAMMIT
    print "******************************"
    print "watch yo damn self, the axis scale '1e7' has been added by hand and may not be correct"
    print "final time =", times[-1]
    print "******************************"
    print n
    ax.text(1.1*n, times[-1], zlim[0]-0.3*(zlim[1] - zlim[0]), '1e8', fontsize=LABELSIZE, zorder=5)
    # scale_txt = ax.yaxis.get_children()[1].get_text()
    # ax.yaxis.get_children()[1].set_position((1, -0.1))
    # ax.yaxis.get_children()[1].set_va('bottom')
    if show:
        fig.tight_layout()
        plt.show()
    return zlim

def plot_degree_surface_v2(degs, times):
    # this method has become overly specialized: it plots many different things, none of which, in fact, is a degree surface
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    n = params['n']
    ci = params['dataInterval']
    indices = np.linspace(1, n, n)
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    #add ones vector
    data_count = 0
    ntimes = times.shape[0]
    max_degs = [np.max(degs[i,:]) for i in range(ntimes)]
    NPLOTTED_PTS = 15
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    FONTSIZE = 20
    LABELSIZE = 16
    fig1 = plt.figure(facecolor='w')
    ax1_ = fig1.add_subplot(111, projection='3d')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    sorted_degs = np.sort(degs, 1)
    #trimmed_degs = sorted_degs
    time_limit = 2*np.power(n, 3)
    interval_times = []
    trimmed_degs = []
    i = 0
    while times[i] <= time_limit:
        interval_times.append(times[i])
        trimmed_degs.append(sorted_degs[i])
        i = i + 1
    interval_times = np.array(interval_times)
    trimmed_degs = np.array(trimmed_degs)
    npoints = trimmed_degs.shape[0]
    NPLOTTED_PTS = 15
    ones = np.ones(NPLOTTED_PTS)
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for v in range(n):
        # ax1_.scatter(100.0*v*ones/n, interval_times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], trimmed_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)],v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v*ones/n, interval_times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], trimmed_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, interval_times[-1], trimmed_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0)
#    plt.show()
    fig2 = plt.figure(facecolor='w')
    ax1_ = fig2.add_subplot(111, projection='3d')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    # sorted_degs = np.log(np.sort(degs, 1)+1)
    sorted_degs = np.sort(degs, 1)
    #trimmed_degs = sorted_degs
    time_limit = 2*np.power(n, 3)
    interval_times = []
    i = 0
    interval_times = np.array(interval_times)
    NPLOTTED_PTS = 15
    ones = np.ones(NPLOTTED_PTS)
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    for v in range(n):
        # ax1_.scatter(100.0*v*ones/n, interval_times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], trimmed_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)],v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v*ones/n, times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], sorted_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, times[-1], sorted_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0)
    #n3
    deg_diff = []
    for i in range(NPLOTTED_PTS):
        deg_diff.append(sorted_degs[i*PLOT_INTERVAL, :] - sorted_degs[(i+1)*PLOT_INTERVAL, :])
    deg_diff = np.array(deg_diff)
    fig13 = plt.figure(facecolor='w')
    ax = fig13.add_subplot(111)
    colornorm = colors.Normalize(vmin=0, vmax=NPLOTTED_PTS-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='RdBu')
    alphaval = 0.1
    for time in range(NPLOTTED_PTS):
        ax.scatter(indices, deg_diff[time], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=alphaval)
    ax.set_xlabel('percentile', fontsize=FONTSIZE)
    ax.set_ylabel('change in degree', fontsize=FONTSIZE)
    ax.set_xlim((0, n))
    ax.set_xticks([i for i in np.linspace(0, n, 11)])
    ax.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    #n2
    deg_diff = []
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for i in range(NPLOTTED_PTS):
        deg_diff.append(sorted_degs[i*PLOT_INTERVAL, :] - sorted_degs[(i+1)*PLOT_INTERVAL, :])
    deg_diff = np.array(deg_diff)
    fig14 = plt.figure(facecolor='w')
    ax = fig14.add_subplot(111)
    for time in range(NPLOTTED_PTS):
        ax.scatter(indices, deg_diff[time], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=alphaval)
    ax.set_xlabel('percentile', fontsize=FONTSIZE)
    ax.set_ylabel('change in degree', fontsize=FONTSIZE)
    ax.set_xlim((0, n))
    ax.set_xticks([i for i in np.linspace(0, n, 11)])
    ax.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    #continue
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    fig6 = plt.figure(facecolor='w')
    # trimmed_degs = sorted_degs
    ax6_ = fig6.add_subplot(111)
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for time in range(NPLOTTED_PTS):
        ax6_.scatter(indices , trimmed_degs[time*PLOT_INTERVAL,:], linewidths=0, c=indices, alpha=1.0*time/NPLOTTED_PTS/5.0)
    ax6_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax6_.set_ylabel('degree', fontsize=FONTSIZE)
    ax6_.set_xlim((0, n))
    ax6_.set_ylim(bottom=0)
    ax6_.set_xticks([i for i in np.linspace(0, n, 11)])
    ax6_.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax6_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    plt.show()
    fig7 = plt.figure(facecolor='w')
    ax6_ = fig7.add_subplot(111)
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    for time in range(NPLOTTED_PTS):
        ax6_.scatter(indices , sorted_degs[time*PLOT_INTERVAL,:], linewidths=0, c=indices, alpha=1.0*time/NPLOTTED_PTS/5.0)
    ax6_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax6_.set_ylabel('degree', fontsize=FONTSIZE)
    ax6_.set_xlim((0, n))
    ax6_.set_ylim(bottom=0)
    ax6_.set_xticks([i for i in np.linspace(0, n, 11)])
    ax6_.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax6_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    plt.show()
#     for v in range(n):
#         ax6_.scatter(100.0*v*ones/n, trimmed_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)],v], linewidths=0, c=colormap.to_rgba(1.0*v), alpha= 1.0*v/n)
#def holdout():
    # fig2 = plt.figure(facecolor='w')
    # ax2_ = fig2.add_subplot(111)
    # ax2_.set_xlabel('simulation step', fontsize=FONTSIZE)
    # ax2_.set_ylabel('max vertex degree', fontsize=FONTSIZE)
    # ax2_.plot(times, max_degs)
    # plt.show()
    #fig 3 is a mess in order to get the colormap and colobars working, requires many of the imports seen at the beginning of the fn
    fig3 = plt.figure(facecolor='w')
    gspec = gs.GridSpec(6,6)
    ax31_ = fig3.add_subplot(gspec[:6,:5])
    maxtime = times[-1]
    ax31_.set_xticks([i*maxtime/10.0 for i in range(11)])
    ax32_ = fig3.add_subplot(gspec[:,5])
    ax31_.set_xlabel('simulation step', fontsize=FONTSIZE)
    ax31_.set_ylabel('vertex degree', fontsize=FONTSIZE)
    ax31_.hold(True)
    artist = []
    for v in range(n):
        ax31_.plot(times, sorted_degs[:,v], c=colormap.to_rgba(1.0*v))
    cb = colorbar.ColorbarBase(ax32_, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax31_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax31_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # cb.set_label('pseudo vertex label')
    fig3.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
    plt.show()
    max_index = 0
    while times[max_index] < 10*np.power(n, 3):
        max_index = max_index + 1
    # fig4 = plt.figure(facecolor='w')
    # ax41_ = fig4.add_subplot(gspec[:6,:5])
    # ax42_ = fig4.add_subplot(gspec[:,5])
    # ax41_.set_xlabel('simulation step', fontsize=FONTSIZE)
    # ax41_.set_ylabel('vertex degree', fontsize=FONTSIZE)
    # ax41_.hold(True)
    # for v in range(n):
    #     ax41_.plot(times[:max_index], sorted_degs[:max_index,v], c=colormap.to_rgba(1.0*v))
    # cb = colorbar.ColorbarBase(ax42_, cmap='jet', norm=colorbarnorm, orientation='vertical')
    # ax41_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    # ax41_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # ax42_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    # ax42_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # # cb.set_label('pseudo vertex label')
    # fig4.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
    # plt.show()
    time_limit = 2*np.power(n, 3)
    interval_times = []
    trimmed_degs = []
    i = 0
    while times[i] <= time_limit:
        interval_times.append(times[i])
        trimmed_degs.append(sorted_degs[i])
        i = i + 1
    interval_times = np.array(interval_times)
    trimmed_degs = np.array(trimmed_degs)
    maxtime = interval_times[-1]
    fig5 = plt.figure(facecolor='w')
    ax51_ = fig5.add_subplot(gspec[:6,:5])
    ax52_ = fig5.add_subplot(gspec[:,5])
    ax51_.set_xlabel('simulation step', fontsize=FONTSIZE)
    ax51_.set_ylabel('vertex degree', fontsize=FONTSIZE)
    ax51_.hold(True)
    for v in range(n):
        ax51_.plot(interval_times, trimmed_degs[:,v], c=colormap.to_rgba(1.0*v))
    ax51_.set_xlim(right=maxtime)
    ax51_.set_xticks([i*maxtime/10.0 for i in range(11)])
    cb = colorbar.ColorbarBase(ax52_, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax51_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax51_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax52_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax52_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # cb.set_label('pseudo vertex label')
    fig5.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
    plt.show()

def plot_time_projection_diff(degs, times):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    gspec = gs.GridSpec(6,6)
    deg_diff = []
    n = params['n']
    ci = params['dataInterval']
    indices = np.linspace(1, n, n)
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    data_count = 0
    ntimes = times.shape[0]
    max_degs = [np.max(degs[i,:]) for i in range(ntimes)]
    NPLOTTED_PTS = 60
    # n3
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    FONTSIZE = 48
    LABELSIZE = 36
    sorted_degs = np.sort(degs)
    time_limit = np.power(n, 3)
    interval_times = []
    trimmed_degs = []
    i = 0
    while times[i] <= time_limit:
        interval_times.append(times[i])
        trimmed_degs.append(sorted_degs[i])
        i = i + 1
    interval_times = np.array(interval_times)
    trimmed_degs = np.array(trimmed_degs)
    npoints = trimmed_degs.shape[0]
    for i in range(NPLOTTED_PTS):
        deg_diff.append(sorted_degs[i*PLOT_INTERVAL, :] - sorted_degs[(i+1)*PLOT_INTERVAL, :])
    deg_diff = np.array(deg_diff)
    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    fig13 = plt.figure(facecolor='w')
    ax = fig13.add_subplot(gspec[:6,:5])
    ax2 = fig13.add_subplot(gspec[:,5])
    colornorm = colors.Normalize(vmin=0, vmax=NPLOTTED_PTS-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    for time in range(NPLOTTED_PTS):
        # ax.scatter(indices, deg_diff[time], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=0.7)
        ax.plot(indices, deg_diff[time], c=colormap.to_rgba(1.0*time))
    colorbarnorm = colors.Normalize(vmin=0, vmax=times[-1])
    cb = colorbar.ColorbarBase(ax2, cmap='jet', norm=colorbarnorm, orientation='vertical', format=formatter)
    # ax2.ticklabel_format(style='sci')
    ax2.yaxis.get_children()[1].set_size(LABELSIZE)
    fig13.text(0.82, 0.94, 'time', fontsize=FONTSIZE-4)
    ax.set_xlabel('percentile', fontsize=FONTSIZE)
    ax.set_ylabel('change in degree', fontsize=FONTSIZE)
    ax.set_xlim((0, n))
    ax.set_xticks([i for i in np.linspace(0, n, 11)])
    ax.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax2.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax2.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    # hard coded magnifications of plot
    # # ax.set_xlim(left=0, right = .68*n)
    # # ax.set_ylim((-500, 1000))
    # ax.set_xlim(left=.68*n, right = n)
    # ax.set_ylim((-4000, 2000))
    # ax.set_xlim((.57*n, n))
    # # ax.set_ylim((-500, 2000))
    # ax.set_ylim((-8000, 4000))
    # n2
    deg_diff = []
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for i in range(NPLOTTED_PTS):
        deg_diff.append(sorted_degs[i*PLOT_INTERVAL, :] - sorted_degs[(i+1)*PLOT_INTERVAL, :])
    deg_diff = np.array(deg_diff)
    fig14 = plt.figure(facecolor='w')
    ax = fig14.add_subplot(gspec[:6,:5])
    ax2 = fig14.add_subplot(gspec[:,5])
    for time in range(NPLOTTED_PTS):
        ax.plot(indices, deg_diff[time], c=colormap.to_rgba(1.0*time))
    colorbarnorm = colors.Normalize(vmin=0, vmax=interval_times[-1])
    cb = colorbar.ColorbarBase(ax2, cmap='jet', norm=colorbarnorm, orientation='vertical', format=formatter)
    # ax2.ticklabel_format(style='sci')
    ax2.yaxis.get_children()[1].set_size(LABELSIZE)
    fig14.text(0.82, 0.94, 'time', fontsize=FONTSIZE-4)
    ax.set_xlabel('percentile', fontsize=FONTSIZE)
    ax.set_ylabel('change in degree', fontsize=FONTSIZE)
    ax.set_xlim((0, n))
    ax.set_xticks([i for i in np.linspace(0, n, 11)])
    ax.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax2.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax2.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # # ax.set_xlim((0, 0.62*n))
    # # ax.set_ylim((-200, 600))
    # ax.set_xlim((.62*n, n))
    # ax.set_ylim((-750, 750))
    # ax.set_xlim((.54*n, n))
    # # ax.set_ylim((-100, 1500))
    # ax.set_ylim((-2100, 1000))
    plt.show()

def plot_time_projection(degs, times, params, fig_title="", cmap='jet', ax_fig=None, ax_cb=None):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs

    fig = plt.figure(facecolor='w')
    gspec = gs.GridSpec(6,6)
    show = False
    if ax_fig is None or ax_cb is None:
        ax_fig = fig.add_subplot(gspec[:6,:5])
        ax_cb = fig.add_subplot(gspec[:,5])
        show = True
    sorted_degs = np.sort(degs, axis=1)
    n = degs.shape[1]
    indices = np.linspace(1, n, n)
    nplotted_pts = 60
    if nplotted_pts > times.shape[0]:
        print "not enough data to plot", nplotted_pts, "points"
        return
    sorted_degs = thin_array(sorted_degs, new_npts=nplotted_pts)
    times.shape = (times.shape[0], 1)
    times = thin_array(times, new_npts=nplotted_pts)
    times.shape = (times.shape[0],)
 
    colornorm = colors.Normalize(vmin=0, vmax=nplotted_pts-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap=cmap)

    for i in range(nplotted_pts):
        ax_fig.plot(indices, sorted_degs[i,:], lw=2, c=colormap.to_rgba(1.0*i), alpha=0.6)
        # ax_fig.plot(indices , sorted_degs[i*PLOT_INTERVAL,:], linewidths=0, c=colormap.to_rgba(1.0*i), alpha=0.3)

    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    colorbarnorm = colors.Normalize(vmin=0, vmax=times[-1])
    cb = colorbar.ColorbarBase(ax_cb, cmap=cmap, norm=colorbarnorm, orientation='vertical', format=formatter)

    SCALESIZE = 24
    ax_cb.yaxis.get_children()[1].set_size(SCALESIZE)
    FONTSIZE = 48
    ax_fig.set_xlabel('vertex', fontsize=FONTSIZE)
    ax_fig.set_ylabel('degree', fontsize=FONTSIZE)
    ax_fig.set_xlim((0, n))
    ax_fig.set_xticklabels([str(int(i)) for i in np.linspace(1, n, 6)])
    # ax_fig.set_xlim(left=)
    # ax_fig.set_yticks([i for i in np.linspace(0, n, 11)])
    # ax_fig.set_yticklabels([str(i) for i in np.linspace(0, 100, 11)/100.0])
    textcolor = ax_fig.get_ymajorticklabels()[0].get_color()
    ax_fig.set_title(fig_title, fontsize=FONTSIZE, color=textcolor)
    ax_cb.text(0.25, 1.04, 'step', fontsize=FONTSIZE-4, color=textcolor)
    LABELSIZE = 36
    ax_fig.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    ax_cb.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    if show:
        plt.show(fig)
    return ax_fig

def plot_degree_projection_superfig(degs, times, ax, colornorm, sort=True, fig_title='', cb_label='degree', FONTSIZE=48):
    import matplotlib.cm as cm
    LABELSIZE = 0.75*FONTSIZE
    LEGENDSIZE = 0.5*FONTSIZE
    ALPHA = 0.8

    n = degs.shape[1]
    npts = times.shape[0]
    npts_kept = int(0.2*npts)
    if sort:
        sorted_degs = np.sort(degs)
    else:
        sorted_degs = degs
    times.shape = (npts, 1)
    new_npts=200
    thinned_sorted_degs = thin_array(sorted_degs, new_npts=new_npts)
    thinned_times = thin_array(times, new_npts=new_npts)
    ones_vect = np.ones(thinned_times.shape[0])

    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    for v in range(n):
        ax.scatter(thinned_times, v*ones_vect, color=colormap.to_rgba(thinned_sorted_degs[:, v]), s=30, lw=0, alpha=ALPHA)

    ax.set_xlabel('step', fontsize=FONTSIZE)
    ax.set_ylabel('vertex', fontsize=FONTSIZE)
    ax.set_xlim((0, thinned_times[-1]))
    ax.set_ylim((0, n))
    ax.set_yticklabels([str(int(i)) for i in np.linspace(1, n, 6)])
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.xaxis.get_children()[1].set_size(LABELSIZE)


def plot_degree_projection(degs, times, sort=True, fig_title='', cb_label='degree', axes=None, FONTSIZE=48, show_colorbar=True):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    import matplotlib.ticker as ticker

    LABELSIZE = 0.75*FONTSIZE
    LEGENDSIZE = 0.5*FONTSIZE
    ALPHA = 0.8

    n = degs.shape[1]
    npts = times.shape[0]
    npts_kept = int(0.2*npts)
    if sort:
        sorted_degs = np.sort(degs)
    else:
        sorted_degs = degs
    times.shape = (npts, 1)
    new_npts=300
    thinned_sorted_degs = thin_array(sorted_degs, new_npts=new_npts)
    thinned_times = thin_array(times, new_npts=new_npts)
    ones_vect = np.ones(thinned_times.shape[0])

    # next_max = np.amax(sorted_degs[:,:-1])
    # log_next_max = np.log(next_max)
    # add one to allow log-taking
    # thinned_sorted_log_degs = np.log(thinned_sorted_degs + 1)

    show = False
    if axes is None:
        fig = plt.figure(facecolor='w')
        gspec = gs.GridSpec(6,6)
        ax = fig.add_subplot(gspec[:,:5], axisbg='white')
        ax_cb = fig.add_subplot(gspec[:, 5])
        show = True
    else:
        ax = axes[0]
        ax_cb = axes[1]

    maxdeg = np.amax(thinned_sorted_degs)
    mindeg = np.amin(thinned_sorted_degs)
    # global colornorm
    colornorm = colors.Normalize(vmin=mindeg, vmax=maxdeg)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    for v in range(n):
        ax.scatter(thinned_times, v*ones_vect, color=colormap.to_rgba(thinned_sorted_degs[:, v]), s=40, lw=0, alpha=ALPHA)

    ax.set_xlabel('Step', fontsize=FONTSIZE)
    ax.set_ylabel('Vertex', fontsize=FONTSIZE)
    ax.set_xlim((0, thinned_times[-1]))
    ax.set_ylim((0, n))
    ax.set_yticklabels([str(int(i)) for i in np.linspace(1, n, 6)])
    ax.tick_params(axis='x', direction='out', top=False)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    ax.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    ax.xaxis.get_children()[1].set_size(LABELSIZE)

    # colorbarnorm = colors.Normalize(vmin=0, vmax=log_next_max)
    # formatting function to prevent long decimal labels
    textcolor = ax.get_ymajorticklabels()[0].get_color()
    ax.set_title(fig_title, fontsize=FONTSIZE, color=textcolor)
    if show_colorbar:
        format_fn = lambda val, pos: "%.3f" % val
        cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colornorm, orientation='vertical', format=ticker.FuncFormatter(format_fn))
        # cb.set_ticks(np.linspace(mindeg, maxdeg, 5))
        ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)
        ax_cb.text(-.2, 1.06, cb_label, fontsize=FONTSIZE-4, color=textcolor)
    if show:
        plt.show()
    return colornorm

def plot_vertex_projection_avg(degs, times):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    n = degs.shape[1]
    degs = np.sort(degs, 1)
    degs = np.average(degs, 2)
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    for i in range(n):
        ax.scatter(times, degs[:, i], c=colormap.to_rgba(1.0*i), lw=0, alpha=0.8)
    ax.set_xlim((0, np.max(times)))
    ax.set_xlabel('step', fontsize=25)
    ax.set_ylabel('vertex degree', fontsize=25)
    ax.tick_params(axis='both', which='both', labelsize=24)
    plt.show()

def plot_vertex_projection(degs, times):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    degs = np.sort(degs)
    n = degs.shape[1]
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    FONTSIZE=48
    LABELSIZE = 0.75*FONTSIZE
    LEGENDSIZE = 0.5*FONTSIZE
    ALPHA = 0.8

    gspec = gs.GridSpec(6,6)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(gspec[:,:5])
    for i in range(n):
        ax.scatter(times, degs[:, i], c=colormap.to_rgba(1.0*i), lw=0, alpha=0.8)

    ax.set_xlim((0, np.max(times)))
    ax.set_xlabel('step', fontsize=FONTSIZE)
    ax.set_ylabel('degree', fontsize=FONTSIZE)
    ax.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    SCALESIZE = 24
    ax.xaxis.get_children()[1].set_size(SCALESIZE)
    # make colorbar
    ax_cb = fig.add_subplot(gspec[:,5])
    cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colornorm, orientation='vertical')
    cb.set_ticks([int(i) for i in np.linspace(0, n-1, 6)])
    cb.set_ticklabels([str(int(i)) for i in np.linspace(1, n, 6)])
    ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    fig.text(0.8, 0.93, 'vertex', fontsize=FONTSIZE-4)
    plt.show()


def plot_vertex_projection_old(degs, times, n3=True):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    gspec = gs.GridSpec(6,6)
    deg_diff = []
    n = params['n']
    ci = params['collection_interval']
    indices = np.linspace(1, n, n)
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    data_count = 0
    ntimes = times.shape[0]
    FONTSIZE = 48
    LABELSIZE = 36
    LEGENDSIZE = 24
    sorted_degs = np.sort(degs)
    # n3
    if n3:
        n3 = np.power(n, 3)
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(gspec[:6,:5])
        maxtime = times[-1]
        n_n3s = float(maxtime/n3)
        n3_indexspacing = int(ntimes/n_n3s)
        n_n3s = int(n_n3s)
        n3_times = times[[i*n3_indexspacing for i in range(1, n_n3s+1)]]
        # ax.set_xticks([i*maxtime/10.0 for i in range(11)])
        ax_cb = fig.add_subplot(gspec[:,5])
        ax.hold(True)
        times.shape = (ntimes, 1)
        thinned_times = thin_array(times, new_npts=ntimes/100)
        thinned_sorted_degs = thin_array(sorted_degs, new_npts=ntimes/100)
        for v in range(n):
            ax.scatter(thinned_times, thinned_sorted_degs[:,v], c=colormap.to_rgba(1.0*v), lw=0, alpha=0.9, s=22.5)
        # ax.axvline(x=n3_times[0], c='0.6', label=r'$n^3$ markers')
        # for i in range(n_n3s):
        #     ax.axvline(x=n3_times[i], c='0.6')
        plot_deviations(ax, n)
        ax.set_xlim((times[0], times[-1]))
        # ax.set_ylim((0, 1.05*np.amax(sorted_degs)))
        ax.set_ylim((0, 2100))
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
        ax.legend(fontsize=LEGENDSIZE, frameon=True, shadow=True)
        cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colorbarnorm, orientation='vertical')
        ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
        ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
        ax_cb.tick_params(axis='both', which='major', labelsize=LABELSIZE)
        ax_cb.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
        # assume ymajorticklabels is not empty
        txtcolor = ax.get_ymajorticklabels()[0].get_color()
        ax.set_xlabel('simulation step', fontsize=FONTSIZE, color=txtcolor)
        ax.set_ylabel('vertex degree', fontsize=FONTSIZE, color=txtcolor)
        ax.set_title(r'$n^3$', fontsize=FONTSIZE, color=txtcolor)
        fig.text(0.78, 0.93, 'percentile', fontsize=FONTSIZE-4, color=txtcolor)
        ax.xaxis.get_children()[1].set_size(LABELSIZE)
        plt.show()
    # n2
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(gspec[:6,:5])
    ax_cb = fig.add_subplot(gspec[:,5])
    n2 = np.power(n, 2)
    time_limit = 20*n2
    interval_times = []
    trimmed_degs = []
    i = 0
    while times[i] <= time_limit:
        interval_times.append(times[i])
        trimmed_degs.append(sorted_degs[i])
        i = i + 1
    interval_times = np.array(interval_times)
    n_n2s = float(interval_times[-1]/n2)
    n2_indexspacing = int(interval_times.shape[0]/n_n2s)
    n_n2s = int(n_n2s)
    n2_times = interval_times[[i*n2_indexspacing for i in range(1, n_n2s+1)]]
    trimmed_degs = np.array(trimmed_degs)
    maxtime = interval_times[-1]
    ax.hold(True)
    for v in range(n):
        ax.scatter(interval_times, trimmed_degs[:,v], c=colormap.to_rgba(1.0*v), lw=0, alpha=0.9, s=25)
    # uncomment for vertical lines marking n2 increments
    # ax.axvline(x=n2_times[0], c='0.6', label=r'$n^2$ markers')
    # for i in range(n_n2s):
    #     ax.axvline(x=n2_times[i], c='0.5')
    plot_deviations(ax, n)
    ax.set_xlim((interval_times[0], interval_times[-1]))
    ax.set_ylim((0, 1.05*np.amax(trimmed_degs)))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    txtcolor = ax.get_ymajorticklabels()[0].get_color()
    ax.legend(fontsize=LEGENDSIZE, frameon=True, shadow=True)
    # ax.set_xticks([i*maxtime/10.0 for i in range(11)])
    cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.set_xlabel('simulation step', fontsize=FONTSIZE, color=txtcolor)
    ax.set_ylabel('vertex degree', fontsize=FONTSIZE, color=txtcolor)
    ax.set_title(r'$n^2$', fontsize=FONTSIZE, color=txtcolor)
    fig.text(0.78, 0.93, 'percentile', fontsize=FONTSIZE-4, color=txtcolor)
    ax.xaxis.get_children()[1].set_size(LABELSIZE)
    plt.show()

def approx_erf(x):
    p = 0.47047
    a1 = 0.3480242
    a2 = -.0958798
    a3 = 0.7478556
    t = lambda x: 1/(1 + p*x)
    tx = t(x)
    return 1 - (a1*tx + a2*np.power(tx, 2) + a3*np.power(tx, 3))*np.exp(-x*x)

def normaldist_cdf(x, mu, sigma):
    return 0.5*(1 - approx_erf((mu - x)/(np.sqrt(2)*sigma)))

def plot_deviations(ax, n):
    dev_plus = lambda T, ndevs: n + ndevs*np.power(2*T/n, 0.5)
    dev_minus = lambda T, ndevs: n - ndevs*np.power(2*T/n, 0.5)
    T0, Tf = ax.get_xlim()
    xvals = np.linspace(T0, Tf, 100)
    devs = np.arange(1, 3)/2.0
    rev_devs = devs[::-1]
    colors = ['k', 'r', 'b', 'w']
    counter = 0
    for ndevs in rev_devs:
        ax.plot(xvals, dev_plus(xvals, ndevs), c=colors[counter], label=str(int(100*0.5*(1 - approx_erf(-ndevs/np.sqrt(2))))/100.0) + ' percentile', lw=3)
        counter = counter + 1
    for ndevs in devs:
        ax.plot(xvals, dev_minus(xvals, ndevs), c=colors[counter], label=str(int(100*0.5*(1 - approx_erf(ndevs/np.sqrt(2))))/100.0) + ' percentile', lw=3)
        counter = counter + 1

def plot_vertex_projection_analytical(degs, times):
    #
    # analytical degree distribution
    #
    # !!!!!!!!!!!!!!! useful only when kappa = 1, rho = 2 !!!!!!!!!!!!!!! 
    #
    # n2
    # time_limit = 2*np.power(n, 3)
    # interval_times = []
    # trimmed_degs = []
    # i = 0
    # sorted_degs = np.sort(degs)
    # while times[i] <= time_limit:
    #     interval_times.append(times[i])
    #     trimmed_degs.append(sorted_degs[i])
    #     i = i + 1
    # interval_times = np.array(interval_times)
    # trimmed_degs = np.array(trimmed_degs)
    # maxtime = interval_times[-1]
    # n3
    ntimes = times.shape[0]
    # reduce times to steady-state values
    # starting_index = int(3.0*ntimes/4.0)
    # plot over entire range for now
    starting_index = 0
    times = times[starting_index:]
    maxtime = times[-1]
    # to thin out number of points plotted
    spacing = 20
    ninittimes = times.shape[0]
    times = np.array([times[i*spacing] for i in range(int(ninittimes/spacing))])
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    gspec = gs.GridSpec(6,6)
    deg_diff = []
    n = params['n']
    ci = params['collection_interval']
    indices = np.linspace(1, n, n)
    FONTSIZE = 48
    LABELSIZE = 36
    # make dataset based on analytical distribution
    kappa = 1
    rho = 2
    poisson_rv = lambda u1, u2: np.random.poisson(4*np.log(1-u1)*np.log(1-u2)/rho)
    nsamples = 10
    sorted_degs = np.zeros(n)
    for k in range(nsamples):
        uniform_indices = np.random.uniform(size=n)
        A_analytical = np.zeros((n,n))
        for i in range(n):
            for j in range(i, n):
                A_analytical[i,j] = poisson_rv(uniform_indices[i], uniform_indices[j])
                A_analytical[j,i] = A_analytical[i,j]
        sorted_degs = sorted_degs + np.sort(np.sum(A_analytical, 0))
    print 'done generating data'
    sorted_degs = sorted_degs/nsamples
    degs = degs[starting_index:,:]
    sorted_degs_fromdata = np.sort(np.array([degs[i*spacing,:] for i in range(int(ninittimes/spacing))]))
    ntimes = times.shape[0]
    poserr = np.max([np.ones(ntimes)*sorted_degs[i] - sorted_degs_fromdata[:,i] for i in range(n)])
    negerr = np.min([np.ones(ntimes)*sorted_degs[i] - sorted_degs_fromdata[:,i] for i in range(n)])
    print poserr, negerr
    errspan = poserr - negerr
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(gspec[:6,:5])
    ax.set_xticks([i*maxtime/10.0 for i in range(11)])
    ax_cb = fig.add_subplot(gspec[:,5])
    ax.set_xlabel('simulation step', fontsize=FONTSIZE)
    ax.set_ylabel('degree difference', fontsize=FONTSIZE)
    ax.hold(True)
    for v in range(n):
        deg_diffs = np.ones(ntimes)*sorted_degs[v]-sorted_degs_fromdata[:,v]
        ax.scatter(times, deg_diffs, c=colormap.to_rgba(1.0*v), lw=0, s=3, alpha=0.3)
    ax.set_xlim((times[0], times[-1]))
    ax.set_ylim((negerr, poserr))
    cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    fig.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
    plt.show()

def plot_densities(densities, times, params):
    n = params['n']
    n_ncubeds = float(times[-1]/np.power(n, 3))
    spacing = int(times.shape[0]/n_ncubeds)
    n_ncubeds = int(n_ncubeds)
    textsize=32
    ticksize=24
    if n_ncubeds > 0:
        ncubed_times = np.array([i*spacing for i in range(1, n_ncubeds+1)]) - 1
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        # n3
        ax.scatter(times, densities, lw=0, alpha=0.7)
        ax.scatter(times[ncubed_times], densities[ncubed_times], c='r', lw=20, marker='_', alpha=0.7, label=r'$n^3$ markers')
        ax.axvline(x=10*np.power(n, 2), c='0.4', label='End of previous graph')
        ax.set_title(r'$n^3$ timescale', fontsize=textsize)
        ax.set_xlabel('Time', fontsize=textsize)
        ax.set_ylabel('Simple edge density', fontsize=textsize)
        ax.legend(fontsize=textsize)
        ax.tick_params(axis='both', which='both', labelsize=ticksize)
        ax.xaxis.get_children()[1].set_size(ticksize)
        ax.set_xlim((times[0], times[-1]))
        ax.set_ylim((0, 1))
        plt.show()
    # n2
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    time_limit = 5*(params['proj_step'] + params['nms'])
    times_short = np.array(filter(lambda t: t <= time_limit, times))
    nshorttimes = times_short.shape[0]
    ax.scatter(times_short, densities[:nshorttimes], lw=0, alpha=0.7, label='simulation')
    plot_densities_analytic(times_short, params, ax, label='analytic')
    ax.set_title(r'$n^2$ timescale', fontsize=textsize)
    ax.set_xlabel('Time', fontsize=textsize)
    ax.set_ylabel('Simple edge density', fontsize=textsize)
    ax.set_xlim((times[0], times_short[-1]))
    ax.set_ylim((0, 1))

    ax.legend(fontsize=textsize)
    ax.tick_params(axis='both', which='both', labelsize=ticksize)
    ax.xaxis.get_children()[1].set_size(ticksize)
    plt.show()

def plot_densities_analytic(times, params, ax, **kwargs):
    n = params['n']
    d = lambda T: 1 - (1 - np.exp(-2*T/(n*n)))*np.exp(-(1 - np.exp(-2*T/(n*n))))
    ax.plot(times, d(times), c='g', lw=3, **kwargs)

def plot_degrees(degs, ax=None):
    if ax is None:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    n = degs.shape[0]
    ax.scatter(np.arange(1, n + 1), -1*np.sort(-degs)/n, lw=0, alpha=0.7)
    ax.plot(np.arange(1, n + 1), -np.log(np.arange(1, n + 1)/float(n)))
    ax.set_ylim((0, np.max(degs[degs != np.max(degs)])/float(n)))
    plt.show()

def comp_selfloops(selfloops, times, params, ax, **kwargs):
    # nms = params['nms']
    # fig = plt.figure(facecolor='w')
    # ax = fig.add_subplot(111)
    # colors = ['b', 'r', 'g', 'k', 'c']
    nruns = len(selfloops)
    n = selfloops[0].shape[0]
    avg_selfloops = np.zeros(n)
    for sl in selfloops:
        avg_selfloops = avg_selfloops + sl
    avg_selfloops = avg_selfloops/nruns
    ax.scatter(times, avg_selfloops, **kwargs)
    #     ax.scatter(times[i], selfloops[i], c=colors[i])
    # ntimes = times[0].shape[0]
    # time = times[0]
    # j = 1
    # for i in range(ntimes-1):
    #     if time[i] < j*nms and time[i+1] > j*nms:
    #         ax.axvline(x=(time[i] + time[i+1])/2.0, c='k')
    #         j = j + 1
    # ax.set_xlim((0, np.max(time)))
    # ax.set_xlabel('step', fontsize=24)
    # ax.set_ylabel('selfloop density', fontsize=24)
    # ax.tick_params(axis='both', which='both', labelsize=24)
    # plt.show()

def compare_deg_recon(pre_recon, post_recon, poly_coeffs):
    n = pre_recon.shape[1]
    nrecons = pre_recon.shape[0]
    ncoeffs = poly_coeffs.shape[1]
    indices = np.arange(n)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    FONTSIZE=48
    LABELSIZE=.75*FONTSIZE
    for i in (2000,):#(500, 1000, 2000):
        ax.scatter(indices, np.sort(pre_recon[i,:]), lw=0, c='b', s=25)
        ax.set_ylabel('degree', fontsize=FONTSIZE)
        ax.set_xlabel('vertex', fontsize=FONTSIZE)
        ax.set_xlim((0,n))
        ax.set_xticklabels([str(int(j)) for j in np.linspace(1, n, 6)])
        ax.tick_params(axis='both', which='both', labelsize=LABELSIZE)
        # ax.scatter(indices, post_recon[i,:], lw=0, c='r')
        recon = np.zeros(n)
        for j in range(ncoeffs):
            recon = indices*recon + np.ones(n)*poly_coeffs[i,ncoeffs-j-1]
        ax.plot(indices, recon, c='g', lw=4)
        eqn_str = "{:1.2e} + ".format(poly_coeffs[i,0]) +  "{:1.2e}".format(poly_coeffs[i,1]) + r'$x$' + " + " +  "{:1.2e}".format(poly_coeffs[i,2]) + r'$x^2$' + " + " +  "{:1.2e}".format(poly_coeffs[i,3]) + r'$x^3$' + " + " +  "{:1.2e}".format(poly_coeffs[i,4]) + r'$x^4$' + " + " +  "{:1.2e}".format(poly_coeffs[i,5]) + r'$x^5$'
        fig.text(0.05, .95, eqn_str, fontsize=0.7*FONTSIZE)
        plt.show()
        plt.savefig('./deg_cpi_data/comparison' + str(i+1) + '.png')
        ax.cla()
# 59.112032033903844, 2.5954064014015832, -0.12128778079321469, 0.0032341112645995402, -3.8690293051141307e-05, 1.6932667805213894e-07

def deg_recon_discrepancy(times, degs, params):
    n = params['n']
    m = n*n/2
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    print degs.shape
    ax.plot(times, np.sum(degs, 1)/2 - m)
    ax.set_xlabel('step')
    ax.set_ylabel('edge discrepancy')
    plt.show()

def align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, ci=None):
    mutual_times = []
    cpi_indices = []
    nocpi_indices = []
    cpi_index = 0
    # manual search
    if ci is None:
        n = cpi_times.shape[0]
        m = nocpi_times.shape[0]
        j = 0
        for i in range(n):
            while j < m and cpi_times[i] > nocpi_times[j]:
                j += 1
            if j < m:
                if cpi_times[i] == nocpi_times[j]:
                    cpi_indices.append(i)
                    nocpi_indices.append(j)
                    mutual_times.append(cpi_times[i])
    else:
        for t in cpi_times:
            if t % ci == 0:
                if t <= nocpi_times[-1]:
                    cpi_indices.append(cpi_index)
                    nocpi_indices.append(t/ci - 1)
                    mutual_times.append(t)
            cpi_index = cpi_index + 1
    return [np.array(mutual_times), cpi_degs[cpi_indices, :], nocpi_degs[nocpi_indices, :]]


def compare_deg_cpi_timeproj(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    import matplotlib.gridspec as gs

    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, ci=nocpi_ci)
    cpi_degs = np.sort(cpi_degs)
    nocpi_degs = np.sort(nocpi_degs)

    # n3

    # uncomment if you want to plot more than one piece of data on the same
    # figure, but this makes things nearly impossible to interpret
    # gspec = gs.GridSpec(6,6)
    # fig = plt.figure(facecolor='w')
    # ax_fig = fig.add_subplot(gspec[:6,:5])
    # ax_cb = fig.add_subplot(gspec[:,5])
    # plot_time_projection(np.abs(cpi_degs - nocpi_degs)/nocpi_degs, mutual_times, cpi_params, fig_title=r'$n^3$', cmap='jet', ax_fig=ax_fig, ax_cb=ax_cb)

    plot_time_projection(cpi_degs, mutual_times, cpi_params, fig_title=r'$n^3$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    plot_time_projection(nocpi_degs, mutual_times, cpi_params, fig_title=r'$n^3$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    # plt.show(fig)

    # n2

    n = cpi_params['n']
    time_limit = 5*(cpi_params['proj_step'] + cpi_params['nms'])
    i = 0
    while mutual_times[i] < time_limit:
        i = i + 1

    # fig = plt.figure(facecolor='w')
    # ax_fig = fig.add_subplot(gspec[:6,:5])
    # ax_cb = fig.add_subplot(gspec[:,5])
    # plot_time_projection(np.abs(cpi_degs[:i, :] - nocpi_degs[:i, :])/nocpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet', ax_fig=ax_fig, ax_cb=ax_cb)

    plot_time_projection(cpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    plot_time_projection(nocpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    # plt.show(fig)

def compare_deg_cpi(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs)#, ci=nocpi_ci)
    cpi_degs = np.sort(cpi_degs)
    nocpi_degs = np.sort(nocpi_degs)

    # implement Mark Brynildsen's suggestion
    plot_total_error(cpi_degs, nocpi_degs, mutual_times)

    # n3
    plot_degree_projection((cpi_degs - nocpi_degs)/nocpi_degs, mutual_times, sort=False, fig_title=r'$n^3$', cb_label='relative error') #, mutual_times, sort=False, fig_title=r'$n^3$', cb_label='relative error') 
    # n2
    n = cpi_params['n']
    time_limit = 5*(cpi_params['proj_step']+cpi_params['nms'])
    i = 0
    while mutual_times[i] < time_limit:
        i = i + 1
    plot_degree_projection((cpi_degs[:i,:] - nocpi_degs[:i,:])/nocpi_degs[:i,:], mutual_times[:i], sort=False, fig_title=r'$n^2$', cb_label='relative error') #, mutual_times[:i], sort=False, fig_title=r'$n^2$', cb_label='relative error') 

def compare_deg_cpi_3d(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    from mpl_toolkits.mplot3d import Axes3D

    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, ci=nocpi_ci)

    # thin arrays, 3D plotting tis expensive
    mutual_times.shape = (mutual_times.shape[0], 1)
    new_npts=100
    cpi_degs = thin_array(cpi_degs, new_npts=new_npts)
    nocpi_degs = thin_array(nocpi_degs, new_npts=new_npts)
    mutual_times = thin_array(mutual_times, new_npts=new_npts)
    mutual_times.shape = (mutual_times.shape[0],)

    n = params['n']
    i = 0
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')
    for t in mutual_times:
        ax.plot(np.linspace(1,n,n), np.ones(n)*t, np.sort(cpi_degs[i,:]), c='r')
        i = i + 1
    ax.plot(np.linspace(1,n,n), np.ones(n)*mutual_times[-1], np.sort(cpi_degs[-1,:]), c='r', label='cpi')
    i = 0
    for t in mutual_times:
        ax.plot(np.linspace(1,n,n), np.ones(n)*t, np.sort(nocpi_degs[i,:]), c='b')
        i = i + 1
    ax.plot(np.linspace(1,n,n), np.ones(n)*mutual_times[-1], np.sort(nocpi_degs[-1,:]), c='b', label='plain')
    fs = 42
    ax.set_xlabel('vertex', fontsize=fs)
    ax.set_ylabel('time', fontsize=fs)
    ax.set_zlabel('degree', fontsize=fs)
    ax.legend(fontsize=fs)
    ax.tick_params(axis='both', which='both', labelsize=fs)
    plt.show()

def newton_deg_evo(deg_seqs):
    import matplotlib.cm as cm
    import matplotlib.colors as colors

    ndeg_seqs = deg_seqs.shape[0]
    colornorm = colors.Normalize(vmin=0, vmax=ndeg_seqs-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    n = deg_seqs.shape[1]
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    for i in range(ndeg_seqs):
        ax.plot(np.arange(n) + 1, np.sort(deg_seqs[i,:]), c=colormap.to_rgba(1.0*i))
    ax.set_xlabel('vertex')
    ax.set_ylabel('degree')
    plt.show()

def newton_resid_evo(resids):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    k = 11
    ax.plot(np.arange(k), resids[:k], zorder=1)
    ax.scatter(np.arange(k), resids[:k], lw=0, c='k', zorder=2)
    ax.set_xlim((0,k-1))
    # ax.plot(np.arange(resids.shape[0]), resids, zorder=1)
    # ax.scatter(np.arange(resids.shape[0]), resids, lw=0, c='k', zorder=2)
    # ax.set_xlim((0,resids.shape[0]-1))
    fs = 48
    ax.tick_params(axis='both', which='both', labelsize=fs)
    ax.set_xlabel('iteration', fontsize=fs)
    ax.set_ylabel(r'$\parallel F(x^{(k)}) \parallel_2$', fontsize=fs)
    ax.tick_params(axis='both', which='both', labelsize=fs)
    plt.show()

def comp_newton(xs, deg_seqs, ax=None):
    if ax is None:
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
    n = xs.shape[1]
    ax.plot(np.arange(n), np.sort(xs[-1,:]), c='b', label='newton solution')
    ax.plot(np.arange(n), np.sort(deg_seqs[-1,:]), c='g', label='direct simulation')
    ax.set_xlabel('vertex', fontsize=24)
    ax.set_ylabel('degree', fontsize=24)
    ax.tick_params(axis='both', which='both', labelsize=24)
    ax.legend(loc=2, fontsize=18)
    plt.show()

def super_comp_newton(xs, deg_seqs, resids):

    import matplotlib.animation as manimation
    FFMpegWriter = manimation.writers['ffmpeg']
    writer = FFMpegWriter(fps=4)

    fig = plt.figure()
    ax_degs = fig.add_subplot(211)
    ax_resid = fig.add_subplot(212)
    n = xs.shape[1]
    npts = xs.shape[0]
    maxdeg = np.amax(deg_seqs)
    mindeg = np.amin(deg_seqs)
    deg_lims = (mindeg - 0.1*(maxdeg-mindeg), maxdeg + 0.1*(maxdeg-mindeg))
    maxresid = np.amax(resids)
    minresid = np.amin(resids)
    resid_lims = (minresid - 0.1*(maxresid-minresid), maxresid + 0.1*(maxresid-minresid))
    
    # anim_fig = plt.figure()
    with writer.saving(fig, "./figs/newton/newton_animation.avi", 100):
        for i in range(npts):
            ax_degs.hold(False)
            ax_degs.plot(np.arange(n), np.sort(xs[i,:]), c='b',  label='newton iteration')
            ax_degs.hold(True)
            ax_degs.plot(np.arange(n), np.sort(deg_seqs[-1,:]), c='g', label='direct simulation')
            ax_degs.set_ylim(deg_lims)
            ax_degs.set_xlabel('vertex', fontsize=24)
            ax_degs.set_ylabel('degree', fontsize=24)
            ax_degs.tick_params(axis='both', which='both', labelsize=24)
            ax_degs.legend(loc=2, fontsize=18)
            ax_resid.hold(False)
            ax_resid.plot(np.arange(1,i+1), resids[:i], zorder=1, c='b')
            ax_resid.hold(True)
            ax_resid.scatter(np.arange(1,i+1), resids[:i], lw=0, c='k', zorder=2)
            ax_resid.set_xlim((1,npts))
            ax_resid.set_ylim(resid_lims)
            ax_resid.set_xlabel('iteration', fontsize=24)
            ax_resid.set_ylabel(r'$\parallel F(x^{(k)}) \parallel_2$', fontsize=24)
            ax_resid.tick_params(axis='both', which='both', labelsize=24)
            fig.tight_layout()
            # some bullshit for ffmpeg
            if i/10 == 0:
                index = str(0) + str(i)
            else:
                index = str(i)
            # wow so good
            if i > 0:
                writer.grab_frame()
            # plt.savefig('./figs/newton/super_comp' + index  + '.png')
        ax_degs.hold(False)
        ax_degs.plot(np.arange(n), np.sort(xs[i,:]), c='b', label='newton iteration')
        ax_degs.hold(True)
        ax_degs.plot(np.arange(n), np.sort(deg_seqs[-1,:]), c='g', label='direct simulation')
        ax_degs.set_ylim(deg_lims)
        ax_degs.set_xlabel('vertex', fontsize=24)
        ax_degs.set_ylabel('degree', fontsize=24)
        ax_degs.tick_params(axis='both', which='both', labelsize=24)
        ax_degs.legend(loc=2, fontsize=18)
        ax_degs.set_ylim(deg_lims)
        ax_resid.plot(np.arange(1, npts+1), resids[:npts], zorder=1, c='b')
        ax_resid.scatter(np.arange(1, npts+1), resids[:npts], lw=0, c='k', zorder=2)
        ax_resid.set_xlim((1,npts))
        ax_resid.set_ylim(resid_lims)
        ax_resid.set_xlabel('iteration', fontsize=24)
        ax_resid.set_ylabel(r'$\parallel F(x^{(k)}) \parallel_2$', fontsize=24)
        ax_resid.tick_params(axis='both', which='both', labelsize=24)
        fig.tight_layout()
        for i in range(npts):
            writer.grab_frame()
            # plt.savefig('./figs/newton/super_comp' + str(npts+i) + '.png')

class FormatAxis:
    def __init__(self, ax, has_zaxis=True):
        self._ax = ax
        self._has_zaxis = has_zaxis

    def format(self, axis, data, format_string, offset=0, nticks=5):
        # dictionary of relevant functions/attributes
        d = {'x':{'ax':self._ax.xaxis, 'set-tick-pos':self._ax.set_xticks, 'set-tick-labels':self._ax.set_xticklabels, 'tick-pos':self._ax.xaxis.get_majorticklocs()},
             'y':{'ax':self._ax.yaxis, 'set-tick-pos':self._ax.set_yticks, 'set-tick-labels':self._ax.set_yticklabels, 'tick-pos':self._ax.yaxis.get_majorticklocs()}}

        if self._has_zaxis:
            d['z'] = {'ax':self._ax.zaxis, 'set-tick-pos':self._ax.set_zticks, 'set-tick-labels':self._ax.set_zticklabels, 'tick-pos':self._ax.get_zticks()}

        # tick positions are constant, regardless of offset
        minval, maxval = np.min(data), np.max(data)
        increment = (maxval - minval)/(nticks-1)
        d[axis]['set-tick-pos']([minval + i*increment for i in range(nticks)])
        # subtract offset from data if using
        if offset != 0:
            if offset < 0:
                offset_str =  '- %1.2f' % abs(offset) # unicode dash u'\u2012'
            else:
                offset_str = '+ %1.2f' % abs(offset)
            # go through the terrible process of figuring out where to put the damn thing
            loc = {'x':0, 'y':0, 'z':0}
            for i,key in enumerate(d.keys()):
                if key is axis:
                    if axis is 'x':
                        loc[key] = np.min(d[key]['tick-pos']) - 0.00*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                    else:
                        loc[key] = np.max(d[key]['tick-pos']) + 0.00*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                else:
                    if key is 'x':
                        loc[key] = np.max(d[key]['tick-pos']) + 0.2*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                    else:
                        loc[key] = np.min(d[key]['tick-pos']) - 0.2*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
            if self._has_zaxis:
                self._ax.text(loc['x'], loc['y'], loc['z'], offset_str, fontsize=12) #maxval-0.05*(maxval-minval)
            else:
                self._ax.text(loc['x'], loc['y'], offset_str, fontsize=12) #maxval-0.05*(maxval-minval)
            data = data - offset
        # set axis tick labels
        minval = np.min(data)
        d[axis]['set-tick-labels']([format_string % (minval + i*increment) for i in range(nticks)])

def pa_dmaps_embedding(eigvals, eigvects, params, t=1, plot_2d=False, plot_3d=True):
    from mpl_toolkits.mplot3d import Axes3D
    ntypes = 2 # params['ntypes']
    n = eigvects.shape[1]
    n_pertype = n/ntypes
    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)[::-1]
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    nvects = eigvals.shape[0]
    output_filename = "pa_embedding_"
    cmaps = ['autumn', 'jet', 'spectral', 'winter', 'summer', 'PrOr', 'RdBu', 'RdYlBu', 'RdYlGn']
    cs = ['b', 'r', 'g', 'y', 'c', 'm', 'k', 'b']
    if plot_2d is False:
        print 'Not plotting 2d embeddings'
    else:
        print 'Plotting 2d embeddings'
        # plot 2d embeddings
        n = eigvects.shape[1]
        for i in [2]: # range(1, nvects):
            for j in [3]: # range(i+1, nvects):
                xvals = np.power(eigvals[i], t)*eigvects[i,:]
                yvals = np.power(eigvals[j], t)*eigvects[j,:]
                fig = plt.figure(facecolor='w')
                ax = fig.add_subplot(111)
                

                xdata = np.empty(ntypes*n_pertype)
                ydata = np.empty(ntypes*n_pertype)

                for k in range(ntypes):
                    xdata[k*n_pertype:(k+1)*n_pertype] = xvals[k*n_pertype:(k+1)*n_pertype]
                    ydata[k*n_pertype:(k+1)*n_pertype] = yvals[k*n_pertype:(k+1)*n_pertype]

                    ax.scatter(xvals[k*n_pertype:(k+1)*n_pertype], yvals[k*n_pertype:(k+1)*n_pertype], c=cs[k], lw=0, alpha=0.3, s=np.arange(1, n_pertype)*1000/(n_pertype), label='trajectory ' + str(k+1))
                    ax.scatter(xvals[(k+1)*n_pertype-1], yvals[(k+1)*n_pertype-1], lw=1, c=cs[k], s=100, alpha=1)
                    # ax.scatter(xvals[k*n_pertype:(k+1)*n_pertype], yvals[k*n_pertype:(k+1)*n_pertype], c=np.arange(n_pertype), lw=0, alpha=0.7, s=np.arange(n_pertype)*np.arange(n_pertype)*100/(n_pertype*n_pertype), cmap=cmaps[k])
                    # ax.scatter(xvals[(k+1)*n_pertype-1], yvals[(k+1)*n_pertype-1], c=0.0, s=100, alpha=0.7, cmap=cmaps[k])
                    ax.hold(True)
                ax.set_xlim((1.05*np.min(xdata), 1.05*np.max(ydata)))
                ax.set_ylim((1.05*np.min(ydata), 1.05*np.max(ydata)))
                ax.set_xlabel('$\\Phi_ ' + str(i) + '$')
                ax.set_ylabel('$\\Phi_ ' + str(j) + '$')
                ax.legend(fontsize=28)
                ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
                ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
                plt.tight_layout()
                # plt.savefig("./figs/embeddings/dmaps/2d/" + output_filename + "eigvects_" + str(i+1) + str(j+1) + ".png")
                plt.show(fig)
                ax.hold(False)
    # plot 3d embeddings
    if plot_3d is False:
        print 'Not plotting 3d embeddings'
    else:
        print 'Plotting 3d embeddings'
        # for i in range(1, nvects):
            # for j in range(i+1, nvects):
                # for k in range(j+1, nvects):
        for i in range(1, 2):
            for j in range(2,3):
                for k in range(3,4):
                    fig = plt.figure(facecolor='w')
                    ax = fig.add_subplot(111, projection='3d')
                    xvals = np.power(eigvals[i], t)*eigvects[i,:]
                    yvals = np.power(eigvals[j], t)*eigvects[j,:]
                    zvals = np.power(eigvals[k], t)*eigvects[k,:]

                    for p in range(ntypes):
                        # ax.scatter(xvals[p*n_pertype:(p+1)*n_pertype], yvals[p*n_pertype:(p+1)*n_pertype], zvals[p*n_pertype:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(n_pertype)*np.arange(n_pertype)*100/(n_pertype*n_pertype))
                        ax.scatter(xvals[p*n_pertype:(p+1)*n_pertype], yvals[p*n_pertype:(p+1)*n_pertype], zvals[p*n_pertype:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(1, n_pertype)*1000/(n_pertype))
                        ax.scatter(xvals[(p+1)*n_pertype-1], yvals[(p+1)*n_pertype-1], zvals[(p+1)*n_pertype-1], lw=50, edgecolor=cs[p], c=cs[p], s=1, alpha=1)
                        # ax.scatter(xvals[p*n_pertype], yvals[p*n_pertype], zvals[p*n_pertype], lw=25, edgecolor=cs[p], c=cs[p], s=1, alpha=1)


                    # ax.hold(False)
                    # ax.scatter(xvals[:n/2], yvals[:n/2], zvals[:n/2], c=np.arange(n/2), lw=0, alpha=1, cmap='Reds')
                    # ax.hold(True)
                    # ax.scatter(xvals[-n/2:], yvals[-n/2:], zvals[-n/2:], c=np.arange(n/2), lw=0, alpha=1, cmap='jet')
                    # ax.scatter(xvals, yvals, zvals, c=np.arange(n), lw=0, alpha=1, cmap='jet')
                    ax.set_xlabel('\n\n\n' + r'$\Phi_ ' + str(i) + '$')
                    ax.set_ylabel('\n\n\n' + r'$\Phi_ ' + str(j) + '$')
                    ax.set_zlabel('\n\n\n' + r'$\Phi_ ' + str(k) + '$')
                    # ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
                    # ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
                    # ax.zaxis.get_major_formatter().set_powerlimits((0, 2))
                    # ax.set_xlim((np.min(xvals), np.max(xvals)))
                    # ax.set_ylim((np.min(yvals), np.max(yvals)))
                    # ax.set_zlim((np.min(zvals), np.max(zvals)))
                    ax.xaxis._axinfo['label']['space_factor'] = 4
                    ax.yaxis._axinfo['label']['space_factor'] = 4
                    ax.zaxis._axinfo['label']['space_factor'] = 4
                    formatter = FormatAxis(ax)
                    formatter.format('x', xvals, '%1.2f', nticks=3)
                    formatter.format('y', yvals, '%1.2f', nticks=3)
                    formatter.format('z', zvals, '%1.2f', nticks=3)
                    # plt.tight_layout()
                    plt.show(fig)
                    # plt.savefig("./figs/embeddings/dmaps/3d/" + output_filename + "eigvects_" + str(i+1) + str(j+1) + str(k+1) + ".png")

                    # probably delete the following to end of fn: is paper-specific 

                    fig = plt.figure(facecolor='w')
                    ax = fig.add_subplot(111, projection='3d')
                    xvals = np.power(eigvals[i], t)*eigvects[i,:]
                    yvals = np.power(eigvals[j], t)*eigvects[j,:]
                    zvals = np.power(eigvals[k], t)*eigvects[k,:]

                    print xvals.shape
                    q = 100
                    xdata = np.empty(2*q)
                    ydata = np.empty(2*q)
                    zdata = np.empty(2*q)
                    for p in range(ntypes):
                        # ax.scatter(xvals[p*n_pertype:(p+1)*n_pertype], yvals[p*n_pertype:(p+1)*n_pertype], zvals[p*n_pertype:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(n_pertype)*np.arange(n_pertype)*100/(n_pertype*n_pertype))
                        xdata[q*p:q*(p+1)] = xvals[(p+1)*n_pertype-q:(p+1)*n_pertype]
                        ydata[q*p:q*(p+1)] = yvals[(p+1)*n_pertype-q:(p+1)*n_pertype]
                        zdata[q*p:q*(p+1)] = zvals[(p+1)*n_pertype-q:(p+1)*n_pertype]
                        ntake = 50
                        avg_finalx = np.average(xvals[(p+1)*n_pertype-ntake:(p+1)*n_pertype])
                        avg_finaly = np.average(yvals[(p+1)*n_pertype-ntake:(p+1)*n_pertype])
                        avg_finalz = np.average(zvals[(p+1)*n_pertype-ntake:(p+1)*n_pertype])
                        

                        ax.scatter(xvals[(p+1)*n_pertype-q:(p+1)*n_pertype], yvals[(p+1)*n_pertype-q:(p+1)*n_pertype], zvals[(p+1)*n_pertype-q:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(1, n_pertype)[-q:]*1000/(n_pertype))
                        # ax.scatter(xvals[(p+1)*n_pertype-1], yvals[(p+1)*n_pertype-1], zvals[(p+1)*n_pertype-1], lw=50, edgecolor=cs[p], c=cs[p], s=1, alpha=1)
                        ax.scatter(avg_finalx, avg_finaly, avg_finalz, lw=75, edgecolor=cs[p], c=cs[p])
                        # ax.scatter(xvals[p*n_pertype], yvals[p*n_pertype], zvals[p*n_pertype], lw=25, edgecolor=cs[p], c=cs[p], s=1, alpha=1)

                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # probably delete the following to match !!! line, is paper-specific
                    

                    # ax.hold(False)
                    # ax.scatter(xvals[:n/2], yvals[:n/2], zvals[:n/2], c=np.arange(n/2), lw=0, alpha=1, cmap='Reds')
                    # ax.hold(True)
                    # ax.scatter(xvals[-n/2:], yvals[-n/2:], zvals[-n/2:], c=np.arange(n/2), lw=0, alpha=1, cmap='jet')
                    # ax.scatter(xvals, yvals, zvals, c=np.arange(n), lw=0, alpha=1, cmap='jet')
                    # formatter = FormatAxis(ax)
                    # formatter.format('x', xdata, '%1.2f', nticks=3)
                    # formatter.format('y', ydata, '%1.3f', nticks=3)
                    # formatter.format('z', zdata, '%1.2f', nticks=3)
                    plt.locator_params(nbins=3)
                    # ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
                    # ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
                    # ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

                    ax.set_xlabel('\n\n\n' + r'$\Phi_ ' + str(i) + '$')
                    ax.set_ylabel('\n\n\n' + r'$\Phi_ ' + str(j) + '$')
                    ax.set_zlabel('\n\n\n' + r'$\Phi_ ' + str(k) + '$')
                    # ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
                    # ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
                    # ax.zaxis.get_major_formatter().set_powerlimits((0, 2))
                    ax.set_xlim((-0.0025, 0.0015))
                    ax.set_zlim((-0.0012, 0.0012))
                    ax.xaxis._axinfo['label']['space_factor'] = 4
                    ax.yaxis._axinfo['label']['space_factor'] = 4
                    ax.zaxis._axinfo['label']['space_factor'] = 4
                    plt.tight_layout()
                    plt.show(fig)
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



def pa_pca_embedding(eigvals, eigvects, orig_data):
    # recall eigenvectors are stored in rows
    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    newbasis_data = np.dot(np.dot(orig_data, np.transpose(eigvects)), eigvects)
    nvects = 10 # eigvals.shape[0] - 1
    output_filename = "pa_embedding_"
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    ax.hold(False)
    n = orig_data.shape[0]
    print eigvects.shape, newbasis_data.shape
    for i in range(nvects):
        for j in range(i+1, nvects):
            ax.scatter(newbasis_data[:, i], newbasis_data[:,j], c=np.arange(n), lw=0, alpha=0.7)
            # ax.set_xlim((np.min(xvals), np.max(xvals)))
            # ax.set_ylim((np.min(yvals), np.max(yvals)))
            ax.set_xlabel('basis vector ' + str(i+1))
            ax.set_ylabel('basis vector ' + str(j+1))
            ax.set_title('pa pca embedding')
            plt.savefig("./figs/embeddings/pca/" + output_filename + "eigvects_" + str(i+1) + str(j+1) + ".png")

def plot_total_error(cpi_degs, nocpi_degs, times):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    scale = np.sum(cpi_degs[0]) # should be constant throughout
    ax.plot(times, np.sum(np.abs(cpi_degs - nocpi_degs), axis=1)/scale)
    ax.set_xlim((0, times[-1]))
    fs = 48
    ax.set_xlabel('Step', fontsize=fs)
    ax.set_ylabel('Relative error in CPI routine', fontsize=fs)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    fig.subplots_adjust(left=0.1, bottom=0.11)
    plt.show()

def comp_final_degs(degs, times, params):
    ntypes = params['ntypes']
    npts = degs.shape[0]
    n = degs.shape[1]
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)
    for i in range(ntypes):
        ax.plot(np.arange(n), degs[npts*(i+1)/ntypes - 1,:])
    plt.show()
            
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    parser.add_argument('-c', '--cmap', nargs=1, default='Paired', choices=[m for m in plt.cm.datad])
    parser.add_argument('--fps', nargs=1, type=int, default=10)
    parser.add_argument('-br', '--bitrate', nargs=1, type=int, default=14400)
    parser.add_argument('-ct', '--container-type', nargs=1, default='.mkv')
    parser.add_argument('--contour', action='store_true', default=False)
    parser.add_argument('--threed', action='store_true', default=False)
    parser.add_argument('--fit', action='store_true', default=False)
    parser.add_argument('--plot-adj-eigvals', action='store_true', default=False)
    parser.add_argument('--plot-adj-leading-eigval', action='store_true', default=False)
    parser.add_argument('--plot-adj-leading-eigvect', action='store_true', default=False)
    parser.add_argument('--plot-lapl-eigvals', action='store_true', default=False)
    parser.add_argument('-svs', '--plot-svs', action='store_true', default=False)
    parser.add_argument('-svEig', '--plot-svsEig', action='store_true', default=False)
    parser.add_argument('-pr', '--plot-reconstruction', action='store_true', default=False)
    parser.add_argument('-ppr', '--parallel-plot-reconstruction', action='store_true', default=False)
    parser.add_argument('-pd', '--plot-degrees', action='store_true', default=False)
    parser.add_argument('-cp', '--compare-projection', '--comp-proj', action='store_true', default=False)
    parser.add_argument('--plot-coeffs', '-coeffs', action='store_true', default=False)
    parser.add_argument('--plot-new-coeffs', action='store_true', default=False)
    parser.add_argument('--plot-name', type=str, default="")
    parser.add_argument('-ds', '--plot-degree-surface', action='store_true', default=False)
    parser.add_argument('-dst', '--ds-time-proj', action='store_true', default=False)
    parser.add_argument('-dsv', '--ds-vertex-proj', action='store_true', default=False)
    parser.add_argument('-dsvavg', '--ds-vertex-proj-avg', action='store_true', default=False)
    parser.add_argument('-dsd', '--ds-degree-proj', action='store_true', default=False)
    parser.add_argument('-dstd', '--ds-time-proj-diff', action='store_true', default=False)
    parser.add_argument('-dsva', '--ds-vertex-proj-analytical', action='store_true', default=False)
    parser.add_argument('--plot-simple-densities', '--plot-sd', action='store_true', default=False)
    parser.add_argument('--plot-selfloop-densities', '--plot-sld', action='store_true', default=False)
    parser.add_argument('--plot-degrees-analytic', '-pda', action='store_true', default=False)
    parser.add_argument('--animate-eigvals', action='store_true', default=False)
    parser.add_argument('--comp-eigvect-recon', action='store_true', default=False)
    parser.add_argument('--comp-deg-recon', action='store_true', default=False)
    parser.add_argument('--comp-selfloops', action='store_true', default=False)
    parser.add_argument('--deg-recon-discrepancy', action='store_true', default=False)
    parser.add_argument('--comp-deg-cpi', action='store_true', default=False)
    parser.add_argument('--comp-deg-cpi-3d', action='store_true', default=False)
    parser.add_argument('--comp-deg-cpi-timeproj', action='store_true', default=False)
    parser.add_argument('--newton', action='store_true', default=False)
    parser.add_argument('--newton-comp', action='store_true', default=False)
    parser.add_argument('--newton-super-comp', action='store_true', default=False)
    parser.add_argument('--dmaps-embeddings', action='store_true', default=False)
    parser.add_argument('--pca-embeddings', action='store_true', default=False)
    parser.add_argument('--super-fig', action='store_true', default=False)
    parser.add_argument('--diffinit-comp', action='store_true', default=False)
    parser.add_argument('--comp-final-degs', action='store_true', default=False)
    args = parser.parse_args()
    # this whole file is a huge piece of
    # the atrocities below won't be noticed
    data_list = []
    params = None
    comp = False
    if args.compare_projection:
        # must have args.inputFiles.size === 2
        fullFile = ''
        cpiFile = ''
        for fileName in args.inputFiles:
            if 'paData_' in fileName:
                fullFile = fileName
            elif 'paDataCPI' in fileName:
                cpiFile = fileName
        paramsFull, dataFull = genData(fullFile)
        paramsCPI, dataCPI = genData(cpiFile)
        compareProjection(dataFull, paramsFull, dataCPI, paramsCPI)
    elif args.comp_final_degs:
        for fileName in args.inputFiles:
            if 'deg' in fileName:
                degs, params = uf.get_data(fileName, header_rows=1)
            elif 'time' in fileName:
                times, params = uf.get_data(fileName, header_rows=1)
        comp_final_degs(degs, times, params)
    elif args.plot_degree_surface or args.ds_time_proj_diff or args.ds_time_proj or args.ds_vertex_proj_analytical or args.ds_vertex_proj or args.ds_degree_proj:
        time_data = []
        deg_data = []
        params = None
        for fileName in args.inputFiles:
            if 'deg' in fileName:
                degs, params = get_data(fileName)
                deg_data.append(degs)
            elif 'time' in fileName:
                time_data, params = get_data(fileName)
        if args.plot_degree_surface:
            # n3
            time_data.shape = (time_data.shape[0], 1)
            # degs = thin_array(deg_data[0], new_npts=50)
            # times = thin_array(time_data, new_npts=50)
            degs = deg_data[0]
            times = time_data
            plot_degree_surface(degs, times, title=r'$n^3$')
            # n2
            n = params['n']
            print "n=", n
            time_limit = np.power(n, 3)/10.0
            i = 0
            while time_data[i] <= time_limit:
                i = i + 1
            degs = thin_array(deg_data[0][:i,:], new_npts=50)
            times = thin_array(time_data[:i], new_npts=50)
            plot_degree_surface(degs, times, title=r'$n^2$')
        if  args.ds_time_proj_diff:
            plot_time_projection_diff(deg_data, time_data)
        if args.ds_time_proj:
            # n3
            plot_time_projection(deg_data[0], time_data, params, r'$n^3$')
            # n2
            n = deg_data[0].shape[0]
            if 'proj_step' in params.keys() and 'nms' in params.keys():
                time_limit = 5*(params['proj_step'] + params['nms'])
            else:
                time_limit = np.power(n, 3)/10.0
            time_limit = 1.2*np.power(n, 2)
            i = 0
            while time_data[i] <= time_limit:
                i = i + 1
            plot_time_projection(deg_data[0][:i, :], time_data[:i], params, r'$n^2$')
        if args.ds_vertex_proj_analytical:
            plot_vertex_projection_analytical(deg_data, time_data)
        if  args.ds_vertex_proj_avg:
            plot_vertex_projection_avg(np.array(deg_data), time_data)
        if  args.ds_vertex_proj:
            plot_vertex_projection(deg_data[0], time_data)
        if args.ds_degree_proj:
            # n3
            plot_degree_projection(deg_data[0], time_data, fig_title=r'$n^3$')
            # n2
            n = params['n']
            time_limit = 5*(params['proj_step'] + params['nms'])
            i = 0
            while time_data[i] <= time_limit:
                i = i + 1
            plot_degree_projection(deg_data[0][:i,:], time_data[:i], fig_title=r'$n^2$')
                
    elif args.plot_simple_densities or args.plot_selfloop_densities:
        density_data = None
        time_data = None
        params = None
        for fileName in args.inputFiles:
            if 'density' in fileName:
                density_data, params = get_data(fileName)
            elif 'time' in fileName:
                time_data, params = get_data(fileName)
        plot_densities(density_data, time_data, params)
    elif args.comp_deg_recon:
        for fileName in args.inputFiles:
            if 'pre' in fileName:
                pre, params = get_data(fileName, header_rows=1)
            elif 'post' in fileName:
                post, params = get_data(fileName, header_rows=1)
            elif 'coeffs' in fileName:
                coeffs, params = get_data(fileName, header_rows=1)
        compare_deg_recon(pre, post, coeffs)
    elif args.deg_recon_discrepancy:
        times = None
        degs = None
        params = None
        for fileName in args.inputFiles:
            if 'times' in fileName:
                times, params = get_data(fileName, header_rows=1)
            elif 'degs' in fileName:
                degs, params = get_data(fileName, header_rows=1)
        deg_recon_discrepancy(times, degs, params)
    elif args.plot_new_coeffs:
        times = None
        coeffs_list = []
        for fileName in args.inputFiles:
            if 'times' in fileName:
                times, params = get_data(fileName, header_rows=1)
            elif 'coeffs' in fileName:
                coeffs, params = get_data(fileName, header_rows=1)
                coeffs_list.append(coeffs)
            # elif 'double' in fileName:
            #     fit, params = get_data(fileName, header_rows=0)
        plot_coeffs(times, coeffs_list, args.plot_name)
        # plot_coeffs_fitting(times, coeffs_list, fit, args.plot_name)
    elif args.comp_deg_cpi or args.comp_deg_cpi_3d or args.comp_deg_cpi_timeproj:
        times_cpi = None
        times_nocpi = None
        degs_cpi = None
        degs_nocpi = None
        cpi_params = None
        nocpi_params = None
        for fileName in args.inputFiles:
            if 'withinit' in fileName:
                if 'times' in fileName:
                    times_cpi, cpi_params = get_data(fileName, header_rows=1)
                elif 'degs' in fileName:
                    degs_cpi, params = get_data(fileName, header_rows=1)
            elif 'noinit' in fileName:
                if 'times' in fileName:
                    times_nocpi, nocpi_params = get_data(fileName, header_rows=1)
                elif 'degs' in fileName:
                    degs_nocpi, params = get_data(fileName, header_rows=1)
            elif 'test' in fileName:
                if 'times' in fileName:
                    times_test, test_params = get_data(fileName, header_rows=1)
                elif 'degs' in fileName:
                    degs_test, params = get_data(fileName, header_rows=1)
        if args.comp_deg_cpi:
            compare_deg_cpi(degs_cpi, times_cpi, cpi_params, degs_nocpi, times_nocpi, nocpi_params)
            compare_deg_cpi(degs_nocpi, times_nocpi, nocpi_params, degs_test, times_test, test_params)
        elif args.comp_deg_cpi_3d:
            compare_deg_cpi_3d(degs_cpi, times_cpi, cpi_params, degs_nocpi, times_nocpi, nocpi_params)
        elif args.comp_deg_cpi_timeproj:
            compare_deg_cpi_timeproj(degs_cpi, times_cpi, cpi_params, degs_nocpi, times_nocpi, nocpi_params)
    elif args.comp_selfloops:
        selfloops_nocpi = []
        times_nocpi = None
        params_nocpi = None
        selfloops_cpi = []
        times_cpi = None
        params_cpi = None
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111)
        for fileName in args.inputFiles:
            if 'noinit' in fileName:
                if 'selfloop_densities' in fileName:
                    sl, params = get_data(fileName, header_rows=1)
                    selfloops_nocpi.append(sl)
                elif 'times' in fileName:
                    times_nocpi, params = get_data(fileName, header_rows=1)
            elif 'withinit' in fileName:
                if 'selfloop_densities' in fileName:
                    sl, params = get_data(fileName, header_rows=1)
                    selfloops_cpi.append(sl)
                elif 'times' in fileName:
                    times_cpi, params = get_data(fileName, header_rows=1)
        comp_selfloops(selfloops_cpi, times_cpi, params_cpi, ax, c='c', label='cpi', lw=0)
        comp_selfloops(selfloops_nocpi, times_nocpi, params_nocpi, ax, c='b', label='no cpi', lw=0)
        ax.legend()
        ax.set_xlim((0, np.max(times_cpi)))
        ax.set_xlabel('step', fontsize=24)
        ax.set_ylabel('selfloop density', fontsize=24)
        ax.tick_params(axis='both', which='both', labelsize=24)
        plt.show()
    elif args.newton:
        for f in args.inputFiles:
            if 'xs' in f:
                xs, g = get_data(f, header_rows=0)
            elif 'resid' in f:
                resids, g = get_data(f, header_rows=0)
        newton_deg_evo(xs)
        newton_resid_evo(resids)
    elif args.newton_comp:
        for f in args.inputFiles:
            if 'xs' in f:
                xs, g = get_data(f, header_rows=0)
            if 'degs' in f:
                degs, g = get_data(f)
        comp_newton(xs, degs)
    elif args.newton_super_comp:
        for f in args.inputFiles:
            if 'xs' in f:
                xs, g = get_data(f, header_rows=0)
            if 'degs' in f:
                degs, g = get_data(f)
            if 'times' in f:
                times, g = get_data(f)
            if 'resids' in f:
                resids, g = get_data(f)
        ax = plot_time_projection(degs, times, g)
        comp_newton(xs, degs, ax)
        super_comp_newton(xs, degs, resids)
    elif args.dmaps_embeddings:
        for f in args.inputFiles:
            if 'eigvects' in f:
                eigvects,  params = get_data(f, header_rows=0)
            if 'eigvals' in f:
                eigvals, params = get_data(f, header_rows=0)
        pa_dmaps_embedding(eigvals, eigvects, params)
    elif args.pca_embeddings:
        for f in args.inputFiles:
            if 'eigvects' in f:
                eigvects, g = get_data(f, header_rows=0)
            elif 'eigvals' in f:
                eigvals, g = get_data(f, header_rows=0)
            elif 'graph_embeddings' in f:
                ges, g = get_data(f, header_rows=0)
        pa_pca_embedding(eigvals, eigvects, ges)
    # elif args.diffinit_comp:
    #     from mpl_toolkits.mplot3d import Axes3D
    #     for f in args.inputFiles:
    #         if 'degs1' in f:
    #             degs1, g = get_data(f)
    #         elif 'degs2' in f:
    #             degs2, g = get_data(f)
    #         elif 'times1' in f:
    #             times1, g = get_data(f)
    #         elif 'times2' in f:
    #             times2, g = get_data(f)
    #     mutual_times, mdegs1, mdegs2 = align_degs(times1, degs1, times2, degs2)
    #     print mutual_times.shape
    #     i = 0
    #     while times2[i] != mutual_times[-1]:
    #         i += 1
    #     fig = plt.figure(facecolor='w')
    #     ax = fig.add_subplot(111, projection='3d')
    #     plot_degree_surface(mdegs1, mutual_times, ax=ax, c='b')
    #     plot_degree_surface(mdegs2, mutual_times, ax=ax, c='r')
    #     plot_degree_surface(np.ones((degs2.shape[0]-i, degs2.shape[1]))*degs1[-1,:], times2[i:], ax=ax, c='b')
    #     plot_degree_surface(degs2[i:,:], times2[i:], ax=ax, c='r')
    #     plt.show(fig)

    elif args.super_fig:
        for f in args.inputFiles:
            # it's too late
            if 'withinit' in f or 'SUPA' in f:
                if 'degs' in f:
                    cpi_degs, g = get_data(f)
                elif 'times' in f:
                    cpi_times, g = get_data(f)
            elif 'noinit' in f:
                if 'degs' in f:
                    nocpi_degs, g = get_data(f)
                elif 'times' in f:
                    nocpi_times, g = get_data(f)
        import matplotlib.gridspec as gs
        import matplotlib.colorbar as colorbar
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure(facecolor='w')
        gspec = gs.GridSpec(36,12)
        gspec.update(wspace=1.5, hspace=30.0)
        fs = 20

        mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs)
        mutual_times.shape = (mutual_times.shape[0], 1)
        # mutual_times = thin_array(mutual_times, new_npts=50)
        # cpi_degs = thin_array(cpi_degs, new_npts=50)
        # nocpi_degs = thin_array(nocpi_degs, new_npts=50)

        # cpi
        ax_3d = fig.add_subplot(gspec[:12, :6], projection='3d', axisbg='w')
        ax_proj_cpi = fig.add_subplot(gspec[:12, 7:-1])
        # ax_proj_cb = fig.add_subplot(gspec[:12, -1])
        cpi_times.shape = (cpi_times.shape[0], 1)
        # cpi_degs_3d = thin_array(cpi_degs, new_npts=50)
        # cpi_times_3d = thin_array(cpi_times, new_npts=50)
        zlim = plot_degree_surface(cpi_degs, mutual_times, ax=ax_3d, FONTSIZE=fs, zlabel=r'$d_{CPI}$')
        cn = plot_degree_projection(cpi_degs, mutual_times, axes=[ax_proj_cpi, None], FONTSIZE=fs, show_colorbar=False)
        # hide axis labels
        ax_proj_cpi.set_xlabel('')
        ax_proj_cpi.tick_params(axis='x', which='both', labelsize=0, direction='out', top=False)
        ax_proj_cpi.xaxis.get_children()[1].set_size(0)
        
        # no cpi
        ax_3d = fig.add_subplot(gspec[12:24, :6], projection='3d', axisbg='w')
        ax_proj_nocpi = fig.add_subplot(gspec[12:24, 7:-1])
        # ax_proj_cb = fig.add_subplot(gspec[12:24, -1])
        # nocpi_times.shape = (nocpi_times.shape[0], 1)
        # nocpi_degs_3d = thin_array(nocpi_degs, new_npts=50)
        # nocpi_times_3d = thin_array(nocpi_times, new_npts=50)
        plot_degree_surface(nocpi_degs, mutual_times, ax=ax_3d, FONTSIZE=fs, zlabel=r'$d_p$', colornorm=cn, zlim=zlim)
        plot_degree_projection_superfig(nocpi_degs, mutual_times, ax_proj_nocpi, cn, FONTSIZE=fs, cb_label='')
        ax_proj_nocpi.set_xlabel('')
        ax_proj_nocpi.tick_params(axis='x', which='both', labelsize=0, direction='out', top=False)
        ax_proj_nocpi.xaxis.get_children()[1].set_size(0)

        LABELSIZE=18
        # one large colorbar

        ax_cb = fig.add_subplot(gspec[1:23, -1])
        format_fn = lambda val, pos: "%i" % val
        cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=cn, orientation='vertical', format=ticker.FuncFormatter(format_fn))
        ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)
        ax_cb.text(-.3, 1.04, 'degree', fontsize=20, color='k')

        # the diff between nocpi and cpi
        ax_3d = fig.add_subplot(gspec[24:, :6], projection='3d', axisbg='w')
        ax_proj_diff = fig.add_subplot(gspec[24:, 7:-1])
        ax_proj_cb = fig.add_subplot(gspec[25:-1, -1])
        plot_degree_surface(np.sort(cpi_degs) - np.sort(nocpi_degs), mutual_times, sort=False, ax=ax_3d, FONTSIZE=fs, zlabel=r'$d_{CPI}-d_p$')
        plot_degree_projection(np.sort(cpi_degs) - np.sort(nocpi_degs), mutual_times, sort=False, axes=[ax_proj_diff, ax_proj_cb], FONTSIZE=fs, cb_label='')
        ax_proj_cb.text(-1.5, 1.1, 'degree difference', fontsize=20, color='k')

        fig.subplots_adjust(left=0.01, right=0.93)
        plt.draw()
        plt.show(fig)
    else:
        for fileName in args.inputFiles:
            params, data = genData(fileName)
            if args.contour:
                animateContour(data, params, args.cmap, args.fps, args.bitrate, args.container_type)
            if args.threed:
                animateRawData(data, params, args.fps, args.bitrate, args.container_type)
            if args.plot_adj_eigvals:
                animateVector(data, params, gGP.getAdjEigVals, args.fps, args.bitrate, args.container_type)
            if args.plot_adj_leading_eigval:
                scalarEvolution(data, params, gGP.getAdjLeadingEigVal)
            if args.plot_adj_leading_eigvect:
                animateVector(data, params, gGP.getAdjLeadingEigVect, args.fps, args.bitrate, args.container_type)
            if args.plot_lapl_eigvals:
                animateVector(data, params, gGP.getLaplEigVals, args.fps, args.bitrate, args.container_type)
            if args.plot_svs:
                animateVector(data, params, gGP.getSVs, args.fps, args.bitrate, args.container_type)
            if args.plot_svsEig:
                animateVector(data, params, gGP.getSVLeadingEigVect, args.fps, args.bitrate, args.container_type)
            if args.plot_reconstruction:
                animateReconstruction(data, params, args.fps, args.bitrate, args.container_type)
            if args.parallel_plot_reconstruction:
                p = Pool(processes=2)
                nData = params['nSteps']/params['dataInterval']
                n = params['n']
                #find max degree to set z limits:
                maxDeg = np.amax(data[:,:])
                newFolder = makeFolder('reconstruction')
                xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
                result = p.map(plotReconstruction, [(xgrid, ygrid, data[i*n:(i+1)*n,:], maxDeg, n, newFolder, genFileName('rawData', params, str(i))) for i in range(nData)])
            if args.plot_degrees:
                makeDegSurface(data, params, gGP.getDegrees)
            if args.plot_coeffs:
                plot_vectors_tc(data, params, args.plot_name)
            if args.fit:
                epsilon = 0.1
                fns = []
                fns.append(lambda x,y: x + y)
                plotFittedData(data, params, fns)
            if 'projData' in fileName:
                data_list.append(data)
                comp = True
            if args.animate_eigvals:
                animate_eigvals(data, params)
            if args.comp_eigvect_recon:
                comp_eigvect_recon(data, params, args.plot_name)
            # if 'eigVectData' in fileName:
                # plotEigVectRecon(data, params)
    if comp:
        compare_recon(data_list, params)
