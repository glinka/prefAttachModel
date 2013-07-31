import getGraphProps as gGP
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
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
    print params
    data = np.genfromtxt(pathToFile, delimiter=",", skip_header=1)
    return params, data

def makeFolder(baseName):
    folderCount = 0
    newFolder = baseName + str(folderCount) + '/'
    while os.path.exists(newFolder):
        folderCount = folderCounter + 1
        newFolder = baseName + str(folderCount) + '/'
    os.mkdir(newFolder)
    return newFolder

def genFileName(baseName, params, uniqueID=''):
    fileName = baseName
    for key in params.keys():
        fileName = fileName + '_' + key + '_' + str(params[key])
    return fileName + uniqueID

def animateContour(data, params, fileName, cmap='Paired', fps=10, bitrate=14400, containerType='.mkv'):
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

def animateRawData(data, params, fileName, fps=10, bitrate=14400, containerType='.mkv'):
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
    maxDeg = np.max([np.max(data[(i)*n:(i+1)*n,:n]) for i in range(nData)])
    for ax in spAxes:
        ax.set_xlim3d(left=0, right=100)
        ax.set_ylim3d(bottom=0, top=100)
        ax.set_zlim3d(bottom=0, top=maxDeg)
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    newFolder = makeFolder('animations')
    for i in range(nData):
        for ax in spAxes:
            ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet')
        plt.draw()
        fileName = genFileName('rawData3d', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

def plotFittedData(data, params, fns):
    from mpl_toolkits.mplot3d import Axes3D
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure()
    fig.hold(True)
    plt.text(n/2,1.05*n,'Fitted data', ha='center', fontsize=16)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(-2.0, 45.0)
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    #don't plot each f(x,y) in wireframe
    stride = 10
    newFolder = makeFolder('fittedPlots')
    for i in range(nData):
        fn, xyScaling = gGP.fitXYFunction(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], fns)
        ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet')
        ax.plot_wireframe(xgrid/xyScaling, ygrid/xyScaling, fn(xgrid, ygrid), rstride=stride, cstride=stride)
        plt.draw()
        fileName = genFileName('fittedPlot', params, str(i))
        plt.savefig(newFolder + fileName + '.png')
        plt.show()
        ax.cla()
    

def animateEigVals(data, params, fileName, fn, fps=10, bitrate=14400, containerType='.mkv'):
    import matplotlib.animation as animation
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.figtext(0.5, 0.92, 'Evolution of adjacency eigenvalues in preferential attachment model', ha='center', fontsize=16)
    plt.xlabel('Eigenvalue index')
    plt.ylabel('Eigvenvalue')
    newFolder = makeFolder('eigVals')
    for i in range(nData):    
        ax.cla()
        ax.plot(np.linspace(1,n,n), fn(data[(i)*n:(i+1)*n,:n], n), marker='o', c=[1,0.5,0.5])
        fileName = genFileName('eigVals', params, str(i))
        plt.savefig(newFolder + fileName + '.png')

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
    args = parser.parse_args()
    for fileName in args.inputFiles:
        params, data = genData(fileName)
        if args.contour:
            animateContour(data, params, fileName, args.cmap, args.fps, args.bitrate, args.container_type)
        if args.threed:
            animateRawData(data, params, fileName, args.fps, args.bitrate, args.container_type)
        if args.plot_adj_eigvals:
            animateEigVals(data, params, fileName, gGP.getAdjEigVals, args.fps, args.bitrate, args.container_type)
        if args.plot_lapl_eigvals:
            animateEigVals(data, params, fileName, gGP.getLaplEigVals, args.fps, args.bitrate, args.container_type)
        #plot fitted data, eventually should be moved into different fn prly
        if args.fit:
            epsilon = 0.1
            fns = []
            fns.append(lambda x,y: x+y)
            fns.append(lambda x,y: 1/(x*y + epsilon))
            plotFittedData(data, params, fns)

