import csv
import numpy as np
import matplotlib.pyplot as plt

import biasWindow


def normalize(angle):

    while angle <= -180.:
        angle += 360.
    while angle > 180:
        angle -= 360.

    return angle


def plot(g, PLOT):

    if type(g) == list:
        g = np.array(g)

    y1 = (g[:, 2]/1.0e6)
    y2 = (g[:, 6]/1.0e6)
    t = g[:, 0]/1.0e6

    w1 = np.zeros_like(t)
    w2 = np.zeros_like(t)

    a1 = np.zeros_like(t)
    a2 = np.zeros_like(t)

    for i in xrange(1, t.shape[0]):
        dt = (t[i]-t[i-1])
        w1[i] = normalize(y1[i]-y1[i-1])/dt
        w2[i] = normalize(y2[i]-y2[i-1])/dt

    w1f = biasWindow.centeredMean(w1)
    w2f = biasWindow.centeredMean(w2)

    w1[0] = next(w1f)
    w2[0] = next(w2f)

    for i in xrange(1, t.shape[0]):
        dt = (t[i]-t[i-1])
        w1[i] = next(w1f)
        w2[i] = next(w2f)
        a1[i] = (w1[i]-w1[i-1])/dt
        a2[i] = (w2[i]-w2[i-1])/dt

    print "Theta: [", np.min(y1), ",", np.max(y1), "]", "[",
    print np.min(y2), ",", np.max(y2), "]"

    print "Omega: [", np.mean(w1), ",", np.std(w1), "]", " [",
    print np.mean(w2), ",", np.std(w2), "]\n"

    print "Alpha: [", np.mean(a1), ",", np.std(a1), "]", " [",
    print np.mean(a2), ",", np.std(a2), "]\n"

    if PLOT == 1:
        plt.close()
        plt.figure(1)
        plt.plot(t, y1)
        plt.title("Yaw 1")
        plt.xlabel("Time (secs)")
        plt.ylabel("Yaw (degrees)")

        plt.figure(2)
        plt.plot(t, y2)
        plt.title("Yaw 2")
        plt.xlabel("Time (secs)")
        plt.ylabel("Yaw (degrees)")

        plt.figure(3)
        plt.plot(t, w1)
        plt.title("Yaw rate 1")
        plt.xlabel("Time (secs)")
        plt.ylabel("Yaw rate (degrees/sec)")

        plt.figure(4)
        plt.plot(t, w2)
        plt.title("Yaw rate 2")
        plt.xlabel("Time (secs)")
        plt.ylabel("Yaw rate (degrees/sec)")

        plt.figure(5)
        plt.plot(t, a1)
        plt.title("Angular acceleration 1")
        plt.xlabel("Time (secs)")
        plt.ylabel("Angular acceleration (degrees/(sec^2))")

        plt.figure(6)
        plt.plot(t, a2)
        plt.title("Angular acceleration 2")
        plt.xlabel("Time (secs)")
        plt.ylabel("Angular acceleration (degrees/(sec^2))")

        plt.show()

    elif PLOT == 2:

        plt.close()

        plt.figure(1)
        plt.hist2d(w1, a1, bins=150, normed=True)
        plt.title("Yaw rate change vs yaw rate(1)")
        plt.xlabel("Yaw rate")
        plt.ylabel("Rate of change of yaw rate")
        plt.colorbar()

        plt.figure(2)
        plt.hist2d(w2, a2, bins=150, normed=True)
        plt.title("Yaw rate change vs yaw rate(2)")
        plt.xlabel("Yaw rate")
        plt.ylabel("Rate of change of yaw rate")
        plt.colorbar()

        plt.show()


def parse(files):

    PLOT = 2

    obs = [list() for i in xrange(len(files))]

    for j in xrange(len(files)):

        print "Gyroscope ", j+1
        print "*"*50

        with open(files[j]) as datafile:
            data = csv.reader(datafile)
            data.next()

            count = 0

            for i in xrange(int(20.91e3)):
                data.next()

            for row in data:
                obs[j].append([float(item) for item in row])

                count += 1
                if count % 2.0e6 == 0:
                    print "Until observation ", count

                    plot(obs[j], PLOT)
                    obs[j] = []

            print "Until observation ", count
            plot(obs[j], PLOT)
            # obs[j] = []

            print "\n"

    return obs

if __name__ == "__main__":
    files = ['../data/pimu_1.csv', '../data/pimu_2.csv',
             '../data/pimu_3.csv']

    #files = ['../data/pimu_2.csv']

    obs = parse(files)
