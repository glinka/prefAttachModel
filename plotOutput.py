from multiprocessing import Pool
import getGraphProps as gGP
import numpy as np
import os
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
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
    makeAnimation(fileName, newFolder, outputFilename="reconstruction")

def makeAnimation(inputFilename, inputFolder, outputFilename="output", fps=10, bitrate=1600, containerType='.mkv'):
    from subprocess import call
    us1 = 1
    us2 = 0
    while us1 > 0:
        us2 = us1
        us1 = inputFilename.find("_", us2)
    inputFiles = inputFilename[:us2+1] + "%d" + inputFilename[-4]
    call(['cd', inputFolder])
    call(["ffmpeg", "-i", inputFiles, "-r", str(fps), "-b", str(bitrate), outputFilename+containerType])

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
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
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
    import matplotlib.animation as animation
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.figtext(0.5, 0.92, 'Evolution of adjacency eigenvalues in preferential attachment model', ha='center', fontsize=16)
    plt.xlabel('Eigenvalue index')
    plt.ylabel('Eigvenvalue')
    newFolder = makeFolder('vector')
    for i in range(nData):    
        ax.cla()
        ax.plot(np.linspace(1,n,n), fn(data[(i)*n:(i+1)*n,:n], n), marker='o', c=[1,0.5,0.5])
        fileName = genFileName('eigVals', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

def plot3dData((x, y, z, zlim, xylim, folder)):
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax1 = fig.add_subplot(211, projection='3d')
    ax1.grid(b=False)
    ax1.set_xlim3d(left=0, right=xylim)
    ax1.set_ylim3d(bottom=0, top=xylim)
    ax1.set_zlim3d(bottom=0, top=zlim)
    ax1.scatter(x, y, z, c=z, cmap='jet')
    plt.savefig(folder + fileName + ".png")
    plt.close(fig)

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
    parser.add_argument('--plot-lapl-eigvals', action='store_true', default=False)
    parser.add_argument('-svs', '--plot-svs', action='store_true', default=False)
    parser.add_argument('-svEig', '--plot-svsEig', action='store_true', default=False)
    parser.add_argument('-pr', '--plot-reconstruction', action='store_true', default=False)
    parser.add_argument('-ppr', '--parallel-plot-reconstruction', action='store_true', default=False)
    args = parser.parse_args()
    for fileName in args.inputFiles:
        params, data = genData(fileName)
        if args.contour:
            animateContour(data, params, args.cmap, args.fps, args.bitrate, args.container_type)
        if args.threed:
            animateRawData(data, params, args.fps, args.bitrate, args.container_type)
        if args.plot_adj_eigvals:
            animateVector(data, params, gGP.getAdjEigVals, args.fps, args.bitrate, args.container_type)
        if args.plot_lapl_eigvals:
            animateVector(data, params, gGP.getLaplEigVals, args.fps, args.bitrate, args.container_type)
        if args.plot_svs:
            animateVector(data, params, gGP.getSVs, args.fps, args.bitrate, args.container_type)
        if args.plot_svsEig:
            animateVector(data, params, gGP.getSVLeadingEigVect, args.fps, args.bitrate, args.container_type)
        if args.plot_reconstruction:
            animateReconstruction(data, params, args.fps, args.bitrate, args.container_type)
        if args.parallel_plot_reconstruction:
            p = Pool(processes=12)
            nData = params['nSteps']/params['dataInterval']
            n = params['n']
            #find max degree to set z limits:
            maxDeg = np.amax(data[:,:])
            newFolder = makeFolder('reconstruction')
            xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
            result = p.map(plot3dData, [(xgrid, ygrid, data[i*n:(i+1)*n,:], maxDeg, n, newFolder, genFileName('rawData', params, str(i))) for i in range(nData)])
        if args.fit:
            epsilon = 0.1
            fns = []
            fns.append(lambda x,y: x + y)
            plotFittedData(data, params, fns)
