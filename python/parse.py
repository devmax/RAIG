import numpy as np
import matplotlib.pyplot as plt
import csv


def separate(files):

    data = [list() for i in xrange(len(files))]

    for j in range(len(files)):
        with open(files[j]) as datafile:
            obs = csv.reader(datafile)
            #obs.next()
            for sublist in obs:
                data[j].append([(item) for item in sublist])

    print "Shape is ", len(data), "x", len(data[0])

    return data

if __name__ == "__main__":
    files = ['gyrolog.txt']
    data = separate(files)

    PLOT = False

    if PLOT:
        for i in xrange(len(files)):
            g = np.array(data[i])

            dyaw = g[:, 2]*0.0001309
            t = g[:, -1]/1000.0
            t = t - t[0]

            yaw = np.zeros_like(dyaw)

            for i in xrange(1, dyaw.shape[0]):
                yaw[i] = yaw[i-1] + (t[i]-t[i-1])*dyaw[i-1]

            yaw = (yaw*180/np.pi) % 360

            plt.figure()
            plt.plot(t, dyaw)

            plt.figure()
            plt.plot(t, yaw)
            plt.show()
