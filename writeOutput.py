import getGraphProps as gGP
import numpy as np
import os

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

def writeVectorData(data, params, fn):
    nData = params['nSteps']/params['dataInterval']
    n = params['n']
    newFolder = makeFolder('vectData')
    fileName = ""
    #find y upper and lower limits
    for i in range(nData):
        fileName = genFileName('vectData', params, str(i)) + ".csv"
        f = open(newFolder + fileName, "w")
        xs = np.linspace(1,n,n)
        ys = fn(data[(i)*n:(i+1)*n,:n])
        for i in range(n):
            f.write(str(xs[i]) + '\t')
            f.write(str(ys[i]) + '\n')
        f.close()

if __name__=='__main__':
    import sys
    params, data = genData(sys.argv[1])
    writeVectorData(data, params, gGP.getAdjLeadingEigVect)
