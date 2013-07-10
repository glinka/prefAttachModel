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
        ims.append((plt.pcolormesh(data[(i)*params['n']:(i+1)*params['n'],1:params['n']+1], cmap=cmap),))
    tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
    metadata = dict(artist='alexander holiday')
    writer = animation.FFMpegWriter(fps=fps, bitrate=bitrate, metadata=metadata)
    us = fileName.find('_')
    animFileName = "paAnim" + fileName[us:-4] + containerType
    tcAnim.save(animFileName, writer=writer)

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    parser.add_argument('-c', '--cmap', nargs=1, default='Paired', choices=[m for m in plt.cm.datad])
    parser.add_argument('--fps', nargs='1', type=int, default=10)
    parser.add_argument('-br', '--bitrate', nargs=1, type=int, default=14400)
    parser.add_argument('-ct', '--container-type', nargs=1, default='.mkv')
    args = parser.parse_args()
    print args
    for csv in args.inputFiles:
        params, data = genData(csv)
        animateTimeCourse(data, params, csv, args.cmap, args.fps, args.bitrate, args.container_type)
