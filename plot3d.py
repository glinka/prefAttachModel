import getGraphProps as gGP
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

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

def plot3d(params, data, fileName):
    fig = plt.figure()
    fig.hold(True)
    plt.figtext(0.5, 0.92, 'Wireframe', fontsize=16)
    spAxes = []
    spAxes.append(fig.add_subplot(221, projection='3d'))
    spAxes.append(fig.add_subplot(222, projection='3d'))
    spAxes.append(fig.add_subplot(223, projection='3d'))
    spAxes.append(fig.add_subplot(224, projection='3d'))
    spAxes[0].view_init(-2.0, 45.0)
    spAxes[1].view_init(-2, 135)
    spAxes[2].view_init(45, 225)
    spAxes[0].set_zlim3d(bottom=0, top=16)
    spAxes[0].set_xlim3d(left=0, right=100)
    spAxes[0].set_ylim3d(bottom=0, top=100)
    spAxes[1].set_zlim3d(bottom=0, top=16)
    spAxes[1].set_xlim3d(left=0, right=100)
    spAxes[1].set_ylim3d(bottom=0, top=100)
    spAxes[2].set_zlim3d(bottom=0, top=16)
    spAxes[2].set_xlim3d(left=0, right=100)
    spAxes[2].set_ylim3d(bottom=0, top=100)
    spAxes[3].set_zlim3d(bottom=0, top=16)
    spAxes[3].set_xlim3d(left=0, right=100)
    spAxes[3].set_ylim3d(bottom=0, top=100)
    i = 50
    xgrid, ygrid = np.meshgrid(np.arange(params['n']),np.arange(params['n']))
    for ax in spAxes:
        ax.scatter(xgrid, ygrid, data[(i)*params['n']:(i+1)*params['n'],:params['n']], c=data[(i)*params['n']:(i+1)*params['n'],:params['n']], cmap='jet', alpha=0.8)
    #fitPlane(params, data, ax)
    us = fileName.find('_')
    #bootleg but functional method for getting filename:
    figFileName = 'pa3d' + fileName[us:-4] + '.png'
    plt.show()
    if not os.path.exists('./temp1/'):
        os.makedirs('./temp1')
    plt.savefig('./temp1/' + figFileName)
    
def fitPlane(params, data, ax=False):
    n = params['n']
    nData = params['nSteps']/params['dataInterval']
    X, Y = np.meshgrid(np.arange(params['n']),np.arange(params['n']))
    a = []
    for k in range(nData):
        Adj = data[k*n:(k+1)*n,:n]
        numerator = 0
        for i in range(n):
            for j in range(n):
                numerator = numerator + (X[i,j] + Y[i,j])*Adj[i,j]
        a.append(numerator/sum(sum((X+Y)**2)))
    def z(x, y):
        return [[a[50]*(x[i,j]+y[i,j]) for j in range(n)] for i in range(n)]
        # Z = np.zeros(x.shape)
        # for i in range(n):
        #     for j in range(n):
        #         Z[i,j] = a[50]*(X[i,j]+Y[i,j])
        # return Z
    if ax != False:
        # threshold = 16
        # XTrim = np.array([X[i,:] for i in range(0, n, threshold)])
        # YTrim = np.array([Y[i,:] for i in range(0, n, threshold)])
        ax.plot_wireframe(X, Y, z(X, Y), colors=[0.7, 1, 0.6], alpha=1)
    return a

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', nargs='+')
    args = parser.parse_args()
    for fileName in args.inputFiles:
        params, data = genData(fileName)
        plot3d(params, data, fileName)
