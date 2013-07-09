import numpy as np
import os
#import matplotlib
#matplotlib.use("Agg")

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

if __name__=="__main__":
    import sys
    import matplotlib.pyplot as plt
    #import matplotlib.animation as animation
    params, data = genData(sys.argv[1])
    fig = plt.figure()
    ims = []
    f = open('pngList.txt', "w")
    for i in range(params['nSteps']/params['dataInterval']):
        fig = plt.figure(i)
        plt.pcolormesh(data[(i)*params['n']:(i+1)*params['n'],1:params['n']+1], cmap='Greys')
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

