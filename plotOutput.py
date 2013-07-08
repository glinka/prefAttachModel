import numpy as np
import os
#import matplotlib
#matplotlib.use("Agg")

def genData(fileName):
    pathToFile = os.path.realpath(fileName)
    f = open(pathToFile, "r")
    paramstr = f.readline()
    params = []
    f.close()
    comma1 = 0
    for i in range(4):
        comma2 = paramstr.find(",")
        params.append(float(paramstr[comma1:comma2]))
        paramstr = paramstr[comma2+1:]
    params.append(float(paramstr))
    for i in range(len(params)):
        if(params[i] % 1 == 0):
            params[i] = int(params[i])
    print params
    data = np.genfromtxt(pathToFile, delimiter=",", skip_header=1)
    print data.shape
    return params, data
#output in key:val to add into dict, loop through first line with while(comma2>0):

if __name__=="__main__":
    import sys
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    params, data = genData(sys.argv[1])
    fig = plt.figure()
    ims = []
    n = params[0]
    m = params[1]
    kappa = params[2]
    nSteps = params[3]
    dataInterval = params[4]
    f = open('pngList.txt', "w")
    for i in range(nSteps/dataInterval):
        fig = plt.figure(i)
        plt.pcolormesh(data[(i)*n:(i+1)*n,1:n+1], cmap='Greys')
        plt.savefig("fig" + str(i) + ".png", format="png")
        f.write("fig" + str(i) + ".png\n")
    f.close()
#        ims.append((plt.pcolormesh(data[(i)*n:(i+1)*n,1:n+1], cmap='Greys'),))
#    tcAnim = animation.ArtistAnimation(fig, ims, interval=10)
#    plt.show()    
     #metadata = dict(artist='alexander holiday')
    #writer = FFMpegWriter(fps=5, metadata=metadata)
#    Writer = animation.MencoderWriter(fps=5, codec='ffmpeg', bitrate=1800, metadata={'artist':'alexander holiday'})
#    tcAnim.save('test.mp4', writer=writer)

