import csv

import matplotlib.pyplot as plt
import numpy as numpy


def separate(files):
    obsx = [[] for i in range(len(files))]
    obsy = [[] for i in range(len(files))]
    obsz = [[] for i in range(len(files))]

    temp = [[] for i in range(len(files))]

    for j in range(len(files)):
        with open(files[j]) as datafile:
            obs = csv.reader(datafile)
            obs.next()
            for sublist in obs:
                obsx[j].append(float(sublist[0])*0.0001309)
                obsy[j].append(float(sublist[1])*0.0001309)
                obsz[j].append(float(sublist[2])*0.0001309)
                temp[j].append(float(sublist[3])*0.0001309)

    observations = numpy.array([numpy.array(zip(*obsx)).T,
                                numpy.array(zip(*obsy)).T,
                                numpy.array(zip(*obsz)).T])
    temp = numpy.array(temp)

    print "Shape is ", observations.shape

    return observations, temp

if __name__ == "__main__":
    files = ['big0.csv', 'big1.csv', 'big2.csv', 'big3.csv']
    [obs, temp] = separate(files)

    plt.show()
