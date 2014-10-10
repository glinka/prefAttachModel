from multiprocessing import Pool
import getGraphProps as gGP
import matplotlib.ticker as ticker
import numpy as np
import os

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

# from mpltools import style
# from mpltools import layout

# style.use('ggplot')

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
    fig = plt.figure()
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
    fig = plt.figure()
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
        fig = plt.figure()
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
    fig = plt.figure()
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

def plot_coeffs(times, coeffs_list, plot_name=""):
    coeffs = np.average(np.array(coeffs_list), 0)
    n = coeffs.shape[0]
    ncoeffs = coeffs.shape[1]
    if plot_name is not "":
        plot_name = plot_name + "_"
    for i in range(ncoeffs):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(times, coeffs[:,i], lw=0)
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
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if plot_name is not "":
        plot_name = plot_name + "_"
    for i in range(nyvects):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for j in range(nprojs):
            ax.scatter(data[j*nsaves_per_proj:(j+1)*nsaves_per_proj,nyvects], data[j*nsaves_per_proj:(j+1)*nsaves_per_proj,i], c=(np.sin(float(i)/nyvects), np.cos(float(1-i)/nyvects), 1-float(i)/nyvects), label="coeff: " + str(i+1), lw=0)
        ax.scatter(np.arange(1, nprojs+1)*(nmicrosteps + proj_step), data[np.arange(nprojs)*nsaves_per_proj + nsaves_per_proj - 1, i], lw=0, s=30, c='r')
        ax.set_xlabel('simulation step')
        ax.set_ylabel('coefficient value')
        ax.set_xlim(left=0)
        plt.savefig("coeffs/" + plot_name + "coeff" + str(i) + ".png")
    #ax.legend(loc=6)

def plot_degree_surface(degs, times):
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
    FONTSIZE = 48
    LABELSIZE = 36
    fig1 = plt.figure()
    ax1_ = fig1.add_subplot(111, projection='3d', axisbg='white')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    sorted_degs = np.sort(degs, 1)
    time_limit = 20*np.power(n, 2)
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
    ntrimmedtimes = interval_times.shape[0]
    NPLOTTED_PTS = 45
    plot_interval = int(ntrimmedtimes/NPLOTTED_PTS)
    ones = np.ones(NPLOTTED_PTS)
    # n2
    for v in range(n):
        ax1_.scatter(100.0*v*ones/n, interval_times[[i*plot_interval for i in range(NPLOTTED_PTS)]], trimmed_degs[[i*plot_interval for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, interval_times[-1], trimmed_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0)
    ax1_.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax1_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax1_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax1_.yaxis.get_children()[1].set_size(LABELSIZE)
    fig1.tight_layout()
    fig2 = plt.figure()
    ax1_ = fig2.add_subplot(111, projection='3d', axisbg='white')
    ax1_.set_xlabel('percentile', fontsize=FONTSIZE)
    ax1_.set_ylabel('step', fontsize=FONTSIZE)
    ax1_.set_zlabel('vertex degree', fontsize=FONTSIZE)
    sorted_degs = np.sort(degs, 1)
    # n3
    plot_interval = int(ntimes/NPLOTTED_PTS)
    for v in range(n):
        ax1_.scatter(100.0*v*ones/n, times[[i*plot_interval for i in range(NPLOTTED_PTS)]], sorted_degs[[i*plot_interval for i in range(NPLOTTED_PTS)], v], linewidths=0, c=colormap.to_rgba(1.0*v))
        ax1_.scatter(100.0*v/n, times[-1], sorted_degs[-1, v], linewidths=0, c=colormap.to_rgba(1.0*v))
    ax1_.set_xlim(left=0, right = 100)
    ax1_.set_ylim(bottom=0)
    ax1_.set_zlim(bottom=0, top=2200)
    ax1_.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax1_.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax1_.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax1_.yaxis.get_children()[1].set_size(LABELSIZE)
    fig2.tight_layout()
    plt.show()

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
    fig1 = plt.figure()
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
    fig2 = plt.figure()
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
    fig13 = plt.figure()
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
    fig14 = plt.figure()
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
    fig6 = plt.figure()
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
    fig7 = plt.figure()
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
    # fig2 = plt.figure()
    # ax2_ = fig2.add_subplot(111)
    # ax2_.set_xlabel('simulation step', fontsize=FONTSIZE)
    # ax2_.set_ylabel('max vertex degree', fontsize=FONTSIZE)
    # ax2_.plot(times, max_degs)
    # plt.show()
    #fig 3 is a mess in order to get the colormap and colobars working, requires many of the imports seen at the beginning of the fn
    fig3 = plt.figure()
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
    # fig4 = plt.figure()
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
    fig13 = plt.figure()
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
    fig14 = plt.figure()
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

def plot_time_projection(degs, times, params, fig_title="", cmap='autumn_r', ax_fig=None, ax_cb=None):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs

    fig = plt.figure()
    gspec = gs.GridSpec(6,6)
    show = False
    if ax_fig is None or ax_cb is None:
        ax_fig = fig.add_subplot(gspec[:6,:5])
        ax_cb = fig.add_subplot(gspec[:,5])
        show = True
    sorted_degs = np.sort(degs)
    n = params['n']
    indices = np.linspace(1, n, n)
    nplotted_pts = 60
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
    # ax_fig.set_xlim(left=)
    # ax_fig.set_yticks([i for i in np.linspace(0, n, 11)])
    # ax_fig.set_yticklabels([str(i) for i in np.linspace(0, 100, 11)/100.0])
    textcolor = ax_fig.get_ymajorticklabels()[0].get_color()
    ax_fig.set_title(fig_title, fontsize=FONTSIZE, color=textcolor)
    ax_cb.text(0.25, 1.04, 'time', fontsize=FONTSIZE-4, color=textcolor)
    LABELSIZE = 36
    ax_fig.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)
    ax_cb.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    if show:
        plt.show(fig)

def plot_degree_projection(degs, times, sort=True, fig_title='', cb_label='degree'):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    import matplotlib.ticker as ticker

    FONTSIZE = 48
    LABELSIZE = 36
    LEGENDSIZE = 24
    ALPHA = 0.8

    n = params['n']
    npts = times.shape[0]
    npts_kept = int(0.2*npts)
    if sort:
        sorted_degs = np.sort(degs)
    else:
        sorted_degs = degs
    times.shape = (npts, 1)
    thinned_sorted_degs = thin_array(sorted_degs, new_npts=100)
    thinned_times = thin_array(times, new_npts=100)
    ones_vect = np.ones(thinned_times.shape[0])

    # next_max = np.amax(sorted_degs[:,:-1])
    # log_next_max = np.log(next_max)
    # add one to allow log-taking
    # thinned_sorted_log_degs = np.log(thinned_sorted_degs + 1)

    fig = plt.figure()
    gspec = gs.GridSpec(6,6)
    ax = fig.add_subplot(gspec[:,:5], axisbg='white')

    maxdeg = np.amax(thinned_sorted_degs)
    mindeg = np.min(thinned_sorted_degs)
    colornorm = colors.Normalize(vmin=mindeg, vmax=maxdeg)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    for v in range(n):
        ax.scatter(thinned_times, v*ones_vect, color=colormap.to_rgba(thinned_sorted_degs[:, v]), s=30, lw=0, alpha=ALPHA)

    ax.set_xlabel('step', fontsize=FONTSIZE)
    ax.set_ylabel('vertex', fontsize=FONTSIZE)
    ax.set_xlim((0, thinned_times[-1]))
    ax.set_ylim((0, n))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
    ax.tick_params(axis='both', which='major', labelsize=LABELSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=LABELSIZE)
    ax.xaxis.get_children()[1].set_size(LABELSIZE)

    ax_cb = fig.add_subplot(gspec[:, 5])
    # colorbarnorm = colors.Normalize(vmin=0, vmax=log_next_max)
    # formatting function to prevent long decimal labels
    format_fn = lambda val, pos: "%.2f" % val
    cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colornorm, orientation='vertical', format=ticker.FuncFormatter(format_fn))
    cb.set_ticks(np.linspace(mindeg, maxdeg, 5))
    ax_cb.tick_params(axis='both', which='both', labelsize=LABELSIZE)

    textcolor = ax.get_ymajorticklabels()[0].get_color()
    ax_cb.text(0.0, 1.04, cb_label, fontsize=FONTSIZE-4, color=textcolor)
    ax.set_title(fig_title, fontsize=FONTSIZE, color=textcolor)
    
    plt.show()

def plot_vertex_projection(degs, times):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    fig = plt.figure()
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
        fig = plt.figure()
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
    fig = plt.figure()
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
    fig = plt.figure()
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
        fig = plt.figure()
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
    fig = plt.figure()
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
        fig = plt.figure()
        ax = fig.add_subplot(111)
    n = degs.shape[0]
    ax.scatter(np.arange(1, n + 1), -1*np.sort(-degs)/n, lw=0, alpha=0.7)
    ax.plot(np.arange(1, n + 1), -np.log(np.arange(1, n + 1)/float(n)))
    ax.set_ylim((0, np.max(degs[degs != np.max(degs)])/float(n)))
    plt.show()

def comp_selfloops(selfloops, times, params, ax, **kwargs):
    # nms = params['nms']
    # fig = plt.figure()
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
    for i in range(nrecons):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(indices, np.sort(pre_recon[i,:]), lw=0, c='b')
        ax.scatter(indices, post_recon[i,:], lw=0, c='r')
        recon = np.zeros(n)
        for j in range(ncoeffs):
            recon = indices*recon + np.ones(n)*poly_coeffs[i,ncoeffs-j-1]
        ax.plot(indices, recon, c='g')
        plt.savefig('./deg_cpi_data/comparison' + str(i+1) + '.png')
# 59.112032033903844, 2.5954064014015832, -0.12128778079321469, 0.0032341112645995402, -3.8690293051141307e-05, 1.6932667805213894e-07

def deg_recon_discrepancy(times, degs, params):
    n = params['n']
    m = n*n/2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print degs.shape
    ax.plot(times, np.sum(degs, 1)/2 - m)
    ax.set_xlabel('step')
    ax.set_ylabel('edge discrepancy')
    plt.show()

def align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, nocpi_ci):
    mutual_times = []
    cpi_indices = []
    nocpi_indices = []
    cpi_index = 0
    for t in cpi_times:
        if t % nocpi_ci == 0:
            if t <= nocpi_times[-1]:
                cpi_indices.append(cpi_index)
                nocpi_indices.append(t/nocpi_ci - 1)
                mutual_times.append(t)
        cpi_index = cpi_index + 1
    return [np.array(mutual_times), cpi_degs[cpi_indices, :], nocpi_degs[nocpi_indices, :]]


def compare_deg_cpi_timeproj(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    import matplotlib.gridspec as gs

    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, nocpi_ci)
    cpi_degs = np.sort(cpi_degs)
    nocpi_degs = np.sort(nocpi_degs)

    # n3

    # uncomment if you want to plot more than one piece of data on the same
    # figure, but this makes things nearly impossible to interpret
    # gspec = gs.GridSpec(6,6)
    # fig = plt.figure()
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

    # fig = plt.figure()
    # ax_fig = fig.add_subplot(gspec[:6,:5])
    # ax_cb = fig.add_subplot(gspec[:,5])
    # plot_time_projection(np.abs(cpi_degs[:i, :] - nocpi_degs[:i, :])/nocpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet', ax_fig=ax_fig, ax_cb=ax_cb)

    plot_time_projection(cpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    plot_time_projection(nocpi_degs[:i, :], mutual_times[:i], cpi_params, fig_title=r'$n^2$', cmap='jet') #, ax_cb=ax_cb, ax_fig=ax_fig)
    # plt.show(fig)

def compare_deg_cpi(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, nocpi_ci)
    cpi_degs = np.sort(cpi_degs)
    nocpi_degs = np.sort(nocpi_degs)
    # n3
    plot_degree_projection((cpi_degs - nocpi_degs)/nocpi_degs, mutual_times, sort=False, fig_title=r'$n^3$', cb_label='relative error')
    # n2
    n = cpi_params['n']
    time_limit = 5*(cpi_params['proj_step']+cpi_params['nms'])
    i = 0
    while mutual_times[i] < time_limit:
        i = i + 1
    plot_degree_projection((cpi_degs[:i,:] - nocpi_degs[:i,:])/nocpi_degs[:i,:], mutual_times[:i], sort=False, fig_title=r'$n^2$', cb_label='relative error')

def compare_deg_cpi_3d(cpi_degs, cpi_times, cpi_params, nocpi_degs, nocpi_times, nocpi_params):
    from mpl_toolkits.mplot3d import Axes3D

    nocpi_ci = nocpi_params['collection_interval']
    mutual_times, cpi_degs, nocpi_degs = align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, nocpi_ci)

    # thin arrays, 3D plotting tis expensive
    mutual_times.shape = (mutual_times.shape[0], 1)
    new_npts=100
    cpi_degs = thin_array(cpi_degs, new_npts=new_npts)
    nocpi_degs = thin_array(nocpi_degs, new_npts=new_npts)
    mutual_times = thin_array(mutual_times, new_npts=new_npts)
    mutual_times.shape = (mutual_times.shape[0],)

    n = params['n']
    i = 0
    fig = plt.figure()
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
            plot_degree_surface(deg_data, time_data)
        if  args.ds_time_proj_diff:
            plot_time_projection_diff(deg_data, time_data)
        if args.ds_time_proj:
            # n3
            plot_time_projection(deg_data, time_data, params, r'$n^3$')
            # n2
            n = params['n']
            time_limit = 5*(params['proj_step'] + params['nms'])
            i = 0
            while time_data[i] <= time_limit:
                i = i + 1
            plot_time_projection(deg_data[:i, :], time_data[:i], params, r'$n^2$')
        if args.ds_vertex_proj_analytical:
            plot_vertex_projection_analytical(deg_data, time_data)
        if  args.ds_vertex_proj:
            plot_vertex_projection(np.array(deg_data), time_data)
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
        pre = None
        post = None
        coeffs = None
        for fileName in args.inputFiles:
            if 'pre' in fileName:
                pre, params = get_data(fileName, header_rows=0)
            elif 'post' in fileName:
                post, params = get_data(fileName, header_rows=0)
            elif 'coeffs' in fileName:
                coeffs, params = get_data(fileName, header_rows=0)
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
        plot_coeffs(times, coeffs_list, args.plot_name)
    elif args.comp_deg_cpi or args.comp_deg_cpi_3d or args.comp_deg_cpi_timeproj:
        times_cpi = None
        times_nocpi = None
        degs_cpi = None
        degs_nocpi = None
        cpi_params = None
        nocpi_params = None
        for fileName in args.inputFiles:
            if 'withinit' in fileName or '08' in fileName:
                if 'times' in fileName:
                    times_cpi, cpi_params = get_data(fileName, header_rows=1)
                elif 'degs' in fileName:
                    degs_cpi, params = get_data(fileName, header_rows=1)
            elif 'noinit' in fileName:
                if 'times' in fileName:
                    times_nocpi, nocpi_params = get_data(fileName, header_rows=1)
                elif 'degs' in fileName:
                    degs_nocpi, params = get_data(fileName, header_rows=1)
        if args.comp_deg_cpi:
            compare_deg_cpi(degs_cpi, times_cpi, cpi_params, degs_nocpi, times_nocpi, nocpi_params)
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
        fig = plt.figure()
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
