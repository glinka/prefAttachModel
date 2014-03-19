from multiprocessing import Pool
import getGraphProps as gGP
import numpy as np
import os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    f = open(path_to_file, "r")
    params_str = f.readline()
    params = get_header_data(params_str)
    f.close()
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    print 'done importing data from: ', filename
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
    fig = plt.figure()
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
    fig = plt.figure()
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
    fig = plt.figure()
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
    print nData
    fig = plt.figure()
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

def plotEigVectRecon(data, params):
    n = params['n']
    nData = data.shape[0]/(2)
    fig = plt.figure()
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
    fig = plt.figure()
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
    fig = plt.figure()
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

def scalarEvolution(data, params, fn):
    nData = params['nSteps']/params['dataInterval']
    stepSize = params['dataInterval']
    n = params['n']
    fig = plt.figure()
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
    fig = plt.figure()
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
    fig = plt.figure()
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
    fig = plt.figure()
    fig2 = plt.figure()
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

def plot_vectors_tc(data, params):
    """ Assumes data is arranged as:

    vector1   vector2   vector3 ... vectorN     times
    data[0]   data[1]   data[2]     data[N-2]   data[N-1]

    and plots vectors against times
    """
    nvects = data.shape[1]
    nyvects = nvects - 1
    print data.shape
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hold(False)
    for i in range(nyvects-1):
        ax.scatter(data[:,nyvects], data[:,i], color=(np.sin(i/nyvects), np.cos(1-i/nyvects), 1-i/nyvects), label="coeff: " + str(i+1), lw=0)
        ax.set_xlabel('simulation step')
        ax.set_ylabel('coefficient value')
        ax.set_xlim(left=0)
        plt.savefig("coeffs/coeff" + str(i) + ".png")
    #ax.legend(loc=6)

def plot_degree_surface(degs, times):
    # this method has become overly specialized: it plots many different things, none of which, in fact, is a degree surface
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    gspec = gs.GridSpec(6,6)
    n = params['n']
    ci = params['dataInterval']
    indices = np.linspace(1, n, n)
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    data_count = 0
    ntimes = times.shape[0]
    max_degs = [np.max(degs[i,:]) for i in range(ntimes)]
    NPLOTTED_PTS = 15
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    FONTSIZE = 48
    LABELSIZE = 36
    fig1 = plt.figure()
    ax1_ = fig1.add_subplot(111, projection='3d')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    sorted_degs = np.sort(degs, 1)
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
    # n2
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for v in range(n):
        ax1_.scatter(100.0*v*ones/n, interval_times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], trimmed_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, interval_times[-1], trimmed_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0)
    ax1_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax1_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    fig2 = plt.figure()
    ax1_ = fig2.add_subplot(111, projection='3d')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    sorted_degs = np.sort(degs, 1)
    time_limit = 2*np.power(n, 3)
    NPLOTTED_PTS = 15
    ones = np.ones(NPLOTTED_PTS)
    # n3
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    for v in range(n):
        ax1_.scatter(100.0*v*ones/n, times[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)]], sorted_degs[[i*PLOT_INTERVAL for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, times[-1], sorted_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0)
    ax1_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax1_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
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
    NPLOTTED_PTS = 15
    # n3
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    FONTSIZE = 48
    LABELSIZE = 36
    sorted_degs = np.sort(degs)
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
    for i in range(NPLOTTED_PTS):
        deg_diff.append(sorted_degs[i*PLOT_INTERVAL, :] - sorted_degs[(i+1)*PLOT_INTERVAL, :])
    deg_diff = np.array(deg_diff)
    fig13 = plt.figure()
    ax = fig13.add_subplot(gspec[:6,:5])
    ax2 = fig13.add_subplot(gspec[:,5])
    colornorm = colors.Normalize(vmin=0, vmax=NPLOTTED_PTS-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='RdBu')
    for time in range(NPLOTTED_PTS):
        ax.scatter(indices, deg_diff[time], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=0.3)
    colorbarnorm = colors.Normalize(vmin=0, vmax=times[-1])
    cb = colorbar.ColorbarBase(ax2, cmap='RdBu', norm=colorbarnorm, orientation='vertical')
    fig13.text(0.8, 0.93, 'time', fontsize=FONTSIZE-4)
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
    fig14 = plt.figure()
    ax = fig14.add_subplot(gspec[:6,:5])
    ax2 = fig14.add_subplot(gspec[:,5])
    for time in range(NPLOTTED_PTS):
        ax.scatter(indices, deg_diff[time], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=0.3)
    colorbarnorm = colors.Normalize(vmin=0, vmax=interval_times[-1])
    cb = colorbar.ColorbarBase(ax2, cmap='RdBu', norm=colorbarnorm, orientation='vertical')
    fig14.text(0.8, 0.93, 'time', fontsize=FONTSIZE-4)
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

def plot_time_projection(degs, times):
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
    NPLOTTED_PTS = 15
    data_count = 0
    ntimes = times.shape[0]
    max_degs = [np.max(degs[i,:]) for i in range(ntimes)]
    FONTSIZE = 48
    LABELSIZE = 36
    colorbarnorm = colors.Normalize(vmin=0, vmax=times[-1])
    colornorm = colors.Normalize(vmin=0, vmax=NPLOTTED_PTS-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='RdBu')
    fig7 = plt.figure()
    ax6_ = fig7.add_subplot(gspec[:6,:5])
    ax62_ = fig7.add_subplot(gspec[:,5])
    # n3
    PLOT_INTERVAL = int(ntimes/NPLOTTED_PTS)
    sorted_degs = np.sort(degs)
    for time in range(NPLOTTED_PTS):
        ax6_.scatter(indices , sorted_degs[time*PLOT_INTERVAL,:], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=0.3)
    cb = colorbar.ColorbarBase(ax62_, cmap='RdBu', norm=colorbarnorm, orientation='vertical')
    fig7.text(0.8, 0.93, 'time', fontsize=FONTSIZE-4)
    ax6_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax6_.set_ylabel('degree', fontsize=FONTSIZE)
    ax6_.set_xlim((0, n))
    ax6_.set_ylim(bottom=0)
    ax6_.set_xticks([i for i in np.linspace(0, n, 11)])
    ax6_.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax62_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax62_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax6_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax6_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # # ax6_.set_xlim((0, .68*n))
    # # ax6_.set_ylim((0, 1400))
    # ax6_.set_xlim((.68*n, n))
    # ax6_.set_ylim((1000, 10000))
    # ax6_.set_xlim((0.57*n, n))
    # # ax6_.set_ylim((0, 2100))
    # ax6_.set_ylim((1000, 21000))
    # n2
    interval_times = []
    trimmed_degs = []
    i = 0
    time_limit = 2*np.power(n, 3)
    while times[i] <= time_limit:
        interval_times.append(times[i])
        trimmed_degs.append(sorted_degs[i])
        i = i + 1
    interval_times = np.array(interval_times)
    trimmed_degs = np.array(trimmed_degs)
    npoints = trimmed_degs.shape[0]
    fig6 = plt.figure()
    ax6_ = fig6.add_subplot(gspec[:6,:5])
    ax62_ = fig6.add_subplot(gspec[:,5])
    PLOT_INTERVAL = int(npoints/NPLOTTED_PTS)
    for time in range(NPLOTTED_PTS):
        ax6_.scatter(indices , trimmed_degs[time*PLOT_INTERVAL,:], linewidths=0, c=colormap.to_rgba(1.0*time), alpha=0.3)
    colorbarnorm = colors.Normalize(vmin=0, vmax=interval_times[-1])
    cb = colorbar.ColorbarBase(ax62_, cmap='RdBu', norm=colorbarnorm, orientation='vertical')
    fig6.text(0.8, 0.93, 'time', fontsize=FONTSIZE-4)
    ax6_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax6_.set_ylabel('degree', fontsize=FONTSIZE)
    ax6_.set_xlim((0, n))
    ax6_.set_ylim(bottom=0)
    ax6_.set_xticks([i for i in np.linspace(0, n, 11)])
    ax6_.set_xticklabels([str(i) for i in np.linspace(0, 100, 11)])
    ax62_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax62_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax6_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax6_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    # # ax6_.set_xlim((0, 0.62*n))
    # # ax6_.set_ylim((0, 1200))
    # ax6_.set_xlim((.62*n, n))
    # ax6_.set_ylim((1000, 5500))
    # ax6_.set_xlim((0.54*n, n))
    # # ax6_.set_ylim((0, 2000))
    # ax6_.set_ylim((1500, 12500))
    # plt.show()
    plt.show()

def plot_vertex_projection(degs, times):
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
    # n3
    fig3 = plt.figure()
    ax31_ = fig3.add_subplot(gspec[:6,:5])
    maxtime = times[-1]
    FONTSIZE = 48
    LABELSIZE = 36
    ax31_.set_xticks([i*maxtime/10.0 for i in range(11)])
    ax32_ = fig3.add_subplot(gspec[:,5])
    ax31_.set_xlabel('simulation step', fontsize=FONTSIZE)
    ax31_.set_ylabel('vertex degree', fontsize=FONTSIZE)
    ax31_.hold(True)
    sorted_degs = np.sort(degs)
    for v in range(n):
        ax31_.plot(times, sorted_degs[:,v], c=colormap.to_rgba(1.0*v))
    cb = colorbar.ColorbarBase(ax32_, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax31_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax31_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    fig3.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
    # n2
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
    fig5 = plt.figure()
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
    starting_index = int(3.0*ntimes/4.0)
    times = times[starting_index:]
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
    FONTSIZE = 48
    LABELSIZE = 36
    # make dataset based on analytical distribution
    kappa = 1
    rho = 2
    poisson_rv = lambda u1, u2: np.random.poisson(4*np.log(1-u1)*np.log(1-u2)/rho)
    nsamples = 50
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
    sorted_degs_fromdata = np.sort(degs[starting_index:,:])
    ntimes = times.shape[0]
    poserr = np.max([np.ones(ntimes)*sorted_degs[i] - sorted_degs_fromdata[:,i] for i in range(n)])
    negerr = np.min([np.ones(ntimes)*sorted_degs[i] - sorted_degs_fromdata[:,i] for i in range(n)])
    print poserr, negerr
    errspan = poserr - negerr
    colornorm = colors.Normalize(vmin=0, vmax=n-1)
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')
    fig3 = plt.figure()
    ax31_ = fig3.add_subplot(gspec[:6,:5])
    maxtime = times[-1]
    ax31_.set_xticks([i*maxtime/10.0 for i in range(11)])
    ax32_ = fig3.add_subplot(gspec[:,5])
    ax31_.set_xlabel('simulation step', fontsize=FONTSIZE)
    ax31_.set_ylabel('degree difference', fontsize=FONTSIZE)
    ax31_.hold(True)
    for v in range(n):
        deg_diffs = np.ones(ntimes)*sorted_degs[v]-sorted_degs_fromdata[:,v]
        ax31_.scatter(times, deg_diffs, c=colormap.to_rgba(1.0*v), lw=0, s=3, alpha=0.3)
    ax31_.set_xlim((times[0], times[-1]))
    ax31_.set_ylim((negerr, poserr))
    cb = colorbar.ColorbarBase(ax32_, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax31_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax31_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax32_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    fig3.text(0.8, 0.93, 'percentile', fontsize=FONTSIZE-4)
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
    parser.add_argument('-coeffs', '--plot-coeffs', action='store_true', default=False)
    parser.add_argument('-ds', '--plot-degree-surface', action='store_true', default=False)
    parser.add_argument('-dstd', '--ds-time-proj-diff', action='store_true', default=False)
    parser.add_argument('-dst', '--ds-time-proj', action='store_true', default=False)
    parser.add_argument('-dsva', '--ds-vertex-proj-analytical', action='store_true', default=False)
    parser.add_argument('-dsv', '--ds-vertex-proj', action='store_true', default=False)
    args = parser.parse_args()
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
    elif args.plot_degree_surface or args.ds_time_proj_diff or args.ds_time_proj or args.ds_vertex_proj_analytical or args.ds_vertex_proj:
        time_data = []
        deg_data = []
        params = []
        for fileName in args.inputFiles:
            if 'deg' in fileName:
                deg_data, params = get_data(fileName)
            elif 'time' in fileName:
                time_data, params = get_data(fileName)
        if args.plot_degree_surface:
            plot_degree_surface(deg_data, time_data)
        if  args.ds_time_proj_diff:
            plot_time_projection_diff(deg_data, time_data)
        if args.ds_time_proj:
            plot_time_projection(deg_data, time_data)
        if args.ds_vertex_proj_analytical:
            plot_vertex_projection_analytical(deg_data, time_data)
        if  args.ds_vertex_proj:
            plot_vertex_projection(deg_data, time_data)
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
                plot_vectors_tc(data, params)
            if args.fit:
                epsilon = 0.1
                fns = []
                fns.append(lambda x,y: x + y)
                plotFittedData(data, params, fns)
            if 'projData' in fileName:
                plotCRecon(data, params)
            # if 'eigVectData' in fileName:
                # plotEigVectRecon(data, params)
