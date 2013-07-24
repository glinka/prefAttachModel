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

def animateContour(data, params, fileName, cmap='Paired', fps=10, bitrate=14400, containerType='.mkv'):
    import matplotlib.animation as animation
    fig = plt.figure()
    plt.text(params['n']/2,1.05*params['n'],'Evolution of adjacency matrix in preferential attachment model', ha='center', fontsize=16)
    ims = []
    for i in range(params['nSteps']/params['dataInterval']):
        ims.append((plt.pcolormesh(data[(i)*params['n']:(i+1)*params['n'],:params['n']], cmap=cmap),))
    tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
    metadata = dict(artist='alexander holiday')
    writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate, metadata=metadata)
    us = fileName.find('_')
    animFileName = "paContourAnim" + fileName[us:-4] + containerType
    tcAnim.save(animFileName, writer=writer)

def animate3d(data, params, fileName, fps=10, bitrate=14400, containerType='.mkv'):
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
    fileCounter = 0
    # create new folder for images
    newFolder = './animations' + str(fileCounter) + '/'
    while (os.path.exists(newFolder)) & (fileCounter < 100):
        fileCounter = fileCounter + 1
        newFolder = './animations' + str(fileCounter) + '/'
    os.mkdir(newFolder)
    # plot all data
    for i in range(nData):
        for ax in spAxes:
            ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet')
        plt.draw()
        plt.savefig(newFolder + 'anim3d_fig' + str(i) + '.png')
      # may be faster to save files in directory, then stitch together with ffmpeg. Not much
      # more work at any rate
    # ims = []
    # for i in range(nData):
    #     ims.append((ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet'),))
    # tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
    # metadata = dict(artist='alexander holiday')
    # writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate, metadata=metadata)
    # us = fileName.find('_')
    # animFileName = "pa3dAnim" + fileName[us:-4] + containerType
    # tcAnim.save(animFileName, writer=writer)

def plotFittedData(data, params, fns):
    from mpl_toolkits.mplot3d import Axes3D
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    fig = plt.figure()
    fig.hold(True)
    plt.text(n/2,1.05*n,'Fitted data', ha='center', fontsize=16)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(-2.0, 45.0)
    #makeFolder(name) fn?
    fileCounter = 0
    newFolder = './fittedPlots' + str(fileCounter) + '/'
    while (os.path.exists(newFolder)) & (fileCounter < 100):
        fileCounter = fileCounter + 1
        newFolder = './fittedPlots' + str(fileCounter) + '/'
    os.mkdir(newFolder)
    xgrid, ygrid = np.meshgrid(np.arange(n),np.arange(n))
    stride = 10
    for i in range(nData):
        print i
        fn = gGP.fitXYFunction(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], fns)
        ax.scatter(xgrid, ygrid, data[(i)*n:(i+1)*n,:n], c=data[(i)*n:(i+1)*n,:n], cmap='jet')
        print fn(xgrid, ygrid)
        ax.plot_wireframe(xgrid, ygrid, fn(xgrid, ygrid), rstride=stride, cstride=stride)
        plt.draw()
        plt.savefig(newFolder + 'fittedPlot' + str(i) + '.png')
        plt.show()
        ax.cla()
    

def animateEigValsTimeCourse(data, params, fileName, fn, fps=10, bitrate=14400, containerType='.mkv'):
    import matplotlib.animation as animation
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.figtext(0.5, 0.92, 'Evolution of adjacency eigenvalues in preferential attachment model', ha='center', fontsize=16)
    plt.xlabel('Eigenvalue index')
    plt.ylabel('Eigvenvalue')
    ims = []
    fps = int(params['nSteps']/params['dataInterval']/5.0)
    for i in range(params['nSteps']/params['dataInterval']):
        ims.append((ax.plot(np.linspace(1,params['n'],params['n']), fn(data[(i)*params['n']:(i+1)*params['n'],:params['n']], params['n']), marker='o', c=[1,0.5,0.5])),)
    tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
    metadata = dict(artist='alexander holiday')
    writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate, metadata=metadata)
    us = fileName.find('_')
    #bootleg but functional method for getting filename:
    space1 = str(fn).find(" ")
    space2 = str(fn).find(" ", space1+1)
    prefix = str(fn)[space1+1:space2]
    animFileName = prefix + fileName[us:-4] + containerType
    tcAnim.save(animFileName, writer=writer)

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
    print args
    for fileName in args.inputFiles:
        params, data = genData(fileName)
        if args.contour:
            animateContour(data, params, fileName, args.cmap, args.fps, args.bitrate, args.container_type)
        if args.threed:
            animate3d(data, params, fileName, args.fps, args.bitrate, args.container_type)
        if args.plot_adj_eigvals:
            animateEigValsTimeCourse(data, params, fileName, gGP.getAdjEigVals, args.fps, args.bitrate, args.container_type)
        if args.plot_lapl_eigvals:
            animateEigValsTimeCourse(data, params, fileName, gGP.getLaplEigVals, args.fps, args.bitrate, args.container_type)
        #plot fitted data, eventually should be moved into different fn prly
        if args.fit:
            epsilon = 0.1
            fns = []
            fns.append(lambda x,y: x+y)
            fns.append(lambda x,y: x*y)
            plotFittedData(data, params, fns)

