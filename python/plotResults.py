import numpy as np
import matplotlib.pyplot as plt
import csv
import pdb

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

    t = data[:, 0]
    gt = data[:, 1]
    Bp = data[:, 2]
    mu = data[:, 3]
    Bp_bt = data[:, 4]

    yaw = np.multiply(t, gt)/131.0
    v_yaw = np.multiply(t, Bp)/131.0
    m_yaw = np.multiply(t, mu)/131.0

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

    plt.figure(3)

    plt.plot(yaw, 'g', label="Ground truth")
    plt.plot(v_yaw, 'r', label="Viterbi Estimate")
    plt.plot(m_yaw, 'b', label="Mean estimate")
    plt.title("Yaw values")
    plt.xlabel("Time")
    plt.ylabel("Yaw")

    plt.legend()

    plt.show()
