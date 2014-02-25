import numpy as np
import matplotlib.pyplot as plt
import csv

files = ['/home/dev/Documents/RAIG/data/results']


def read():

    data = list()

    with open(files[0]) as datafile:
        obs = csv.reader(datafile)
        for sublist in obs:
            data.append([float(item) for item in sublist])
        data = np.array(data)

    print "Shape is ", data[0].shape

    return data

if __name__ == '__main__':

    data = read()

    gt = data[:, 0]
    Bp = data[:, 1]
    mu = data[:, 2]
    Bp_bt = data[:, 3]

    plt.figure(1)

    plt.subplot(211)
    plt.plot(gt, 'g', label="Ground truth")
    plt.plot(Bp, 'r', label="Viterbi Estimate")
    plt.plot(mu, 'b', label="Mean estimate")
    plt.title("Instantaneous Estimates")
    plt.xlabel("Time")
    plt.ylabel("Rate")

    plt.legend()

    plt.subplot(212)
    plt.plot(gt, 'g', label="Ground truth")
    plt.plot(Bp_bt, 'r', label="Backtrack Estimate")
    plt.plot(mu, 'b', label="Mean estimate")
    plt.xlabel("Time")
    plt.ylabel("Rate")
    plt.title("Backtracking estimates")

    plt.legend()

    plt.figure(2)

    plt.subplot(211)
    plt.plot(gt-Bp, 'r', label="Viterbi error")
    plt.plot(gt-mu, 'b', label="Mean error")
    plt.xlabel("Time")
    plt.ylabel("Error")

    plt.legend()

    plt.subplot(212)
    plt.plot(gt-Bp_bt, 'r', label="Backtrack error")
    plt.plot(gt-mu, 'b', label="Mean error")
    plt.xlabel("Time")
    plt.ylabel("Error")

    plt.legend()

    plt.show()
