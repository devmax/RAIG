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
            data[j][:, 0] /= 1.0e6
            data[j][:, 1] /= 131.0

    print "Shape is ", len(data), "x", data[0].shape

    return data

if __name__ == "__main__":
    files = ['../data/bias0.txt', '../data/bias1.txt',
             '../data/bias2.txt']
    data = separate(files)

    PLOT = False

    if PLOT:
        for i in xrange(len(files)):
            g = np.array(data[i])

            dyaw = g[:, 1]
            t = g[:, 0]/1.0e6

            yaw = np.zeros_like(dyaw)

            for i in xrange(1, dyaw.shape[0]):
                yaw[i] = yaw[i-1] + (t[i]-t[i-1])*dyaw[i-1]

            yaw = (yaw*180/np.pi) % 360

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
