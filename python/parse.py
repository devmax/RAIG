import numpy as np
import matplotlib.pyplot as plt
import csv


prefix = '/home/dev/Documents/RAIG/data/'
f = ['big0ab', 'big1ab', 'big2ab',
     'big0ac', 'big1ac', 'big2ac',
     'big0ad', 'big1ad', 'big2ad',
     'big0ae', 'big1ae', 'big2ae',
     'big0af', 'big1af', 'big2af',
     'big0ag', 'big1ag', 'big2ag',
     'big0ah', 'big1ah', 'big2ah', ]


def separate(files=f, maxCount=None, maxFiles=9):

    data = [list() for i in xrange(min(len(files), len(files[:maxFiles])))]

    for j in range(len(data)):
        count = 0
        with open(prefix+files[j]) as datafile:
            obs = csv.reader(datafile)
            #obs.next()
            for sublist in obs:
                count += 1
                data[j].append([float(item) for item in sublist])
                if maxCount is not None and count >= maxCount:
                    break
            data[j] = np.array(data[j])

    print "Shape is ", len(data), "x", data[0].shape

    return data

if __name__ == "__main__":

    files = ['xaa', 'xab', 'xac', 'xad', 'xae', 'xaf']

    data = separate(files=files)

    PLOT = True

    if PLOT:
        for i in xrange(len(files)):
            g = np.array(data[i])

            dyaw = g[:1.0e6, 1]
            t = g[:1.0e6, 0]

            plt.figure()
            plt.title("Time vs Rate")
            plt.scatter(t, dyaw)
            #plt.plot(dyaw)

            plt.show()
