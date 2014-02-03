import numpy as np
import matplotlib.pyplot as plt
import csv


def separate(files):

    data = [list() for i in xrange(len(files))]

    for j in range(len(files)):
        count = 0
        with open(files[j]) as datafile:
            obs = csv.reader(datafile)
            #obs.next()
            for sublist in obs:
                count += 1
                data[j].append([float(item) for item in sublist])
                #if count > 1.0e6:
                #    break
            data[j] = np.array(data[j])

    print "Shape is ", len(data), "x", data[0].shape

    return data

if __name__ == "__main__":
    files = ['../data/big2_regaa']

    data = separate(files)

    PLOT = True

    if PLOT:
        for i in xrange(len(files)):
            g = np.array(data[i])

            dyaw = g[:1.0e6, 1]
            t = g[:1.0e6, 0]

            #plt.figure()
            #plt.title("Time vs Rate")
            #plt.plot(t, dyaw)

            plt.figure()
            plt.title("Time")
            plt.plot(t)

            plt.figure()
            plt.title("Rate")
            plt.plot(dyaw)

            plt.show()
