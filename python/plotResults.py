import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb

files = ['/home/dev/Documents/RAIG/data/results']


def read():

    data = list()

    s = 0
    inited = False

    with open(files[0]) as datafile:
        obs = csv.reader(datafile)
        for sublist in obs:
            if inited is False:
                s = len(sublist)
                inited = True
            else:
                if(len(sublist) == s):
                    data.append(([float(item) for item in sublist]))
        data = np.array(data)

    print "Shape is ", data[0].shape

    return data

if __name__ == '__main__':

    data = read()

    dt = data[:, 0]
    gt = data[:, 1]
    v = data[:, 2]
    mu = data[:, 3]
    yaw = data[:, 4]
    v_yaw = data[:, 5]
    mu_yaw = data[:, 6]

    plt.figure(1)

    plt.subplot(211)
    plt.plot(gt, 'g', label="Ground truth")
    plt.plot(v, 'r', label="Viterbi Estimate")
    plt.plot(mu, 'b', label="Mean estimate")
    plt.title("Instantaneous Estimates")
    plt.xlabel("Time")
    plt.ylabel("Rate")

    plt.legend()

    plt.subplot(212)
    plt.plot(abs(gt-v), 'r', label="Viterbi error")
    plt.plot(abs(gt-mu), 'b', label="Mean error")
    plt.xlabel("Time")
    plt.ylabel("Error")

    plt.legend()

    plt.figure(2)

    plt.plot(yaw, 'g', label="Ground truth")
    plt.plot(v_yaw, 'r', label="Viterbi Estimate")
    plt.plot(mu_yaw, 'b', label="Mean estimate")
    plt.title("Yaw values")
    plt.xlabel("Time")
    plt.ylabel("Yaw")

    plt.legend()

    plt.show()
