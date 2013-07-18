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

def animateTimeCourse(data, params, fileName, cmap='Paired', fps=10, bitrate=14400, containerType='.mkv'):
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
    animFileName = "paAnim" + fileName[us:-4] + containerType
    tcAnim.save(animFileName, writer=writer)

def animateEigValsTimeCourse(data, params, fileName, fn):
    import matplotlib.animation as animation
    import ctypes
    fig = plt.figure()
    plt.text(params['n']/2,1.05*params['n'],'Evolution of adjacency eigenvalues in preferential attachment model', ha='center', fontsize=16)
    ims = []
    for i in range(params['nSteps']/params['dataInterval']):
        ptrs = fn(tuple(data[(i)*params['n']:(i+1)*params['n'],:params['n']]))
        doubleList = []
        for j in range(params['n']):
            doubleList.append(cgp.getSnglPtrVal(ptrs, j))
            print cgp.getSnglPtrVal(ptrs, j)
        ims.append((plt.scatter(np.linspace(1,params['n'],params['n']), doubleList[:])))
    tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
    metadata = dict(artist='alexander holiday')
    writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate, metadata=metadata)
    us = fileName.find('_')
    animFileName = "adjEigValsAnim" + fileName[us:-4] + containerType
    tcAnim.save(animFileName, writer=writer)

if __name__=="__main__":
    import argparse
    import calcGraphProps as cgp
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    parser.add_argument('-c', '--cmap', nargs=1, default='Paired', choices=[m for m in plt.cm.datad])
    parser.add_argument('--fps', nargs='1', type=int, default=10)
    parser.add_argument('-br', '--bitrate', nargs=1, type=int, default=14400)
    parser.add_argument('-ct', '--container-type', nargs=1, default='.mkv')
    parser.add_argument('--animate', nargs=1, type=bool, default=False)
    parser.add_argument('--plot-adj-eigvals', nargs=1, type=bool, default=False)
    parser.add_argument('--plot-lapl-eigvals', nargs=1, type=bool, default=False)
    args = parser.parse_args()
    print args
    for fileName in args.inputFiles:
        params, data = genData(fileName)
        if args.animate:
            animateTimeCourse(data, params, fileName, args.cmap, args.fps, args.bitrate, args.container_type)
        if args.plot_adj_eigvals:
            animateEigValsTimeCourse(data, params, fileName, cgp.calcGraphProps_getAdjEigVals)
        if args.plot_lapl_eigvals:
            animateEigValsTimeCourse(data, params, fileName, cgp.calcGraphProps_getLaplEigVals)
            
