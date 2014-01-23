import csv
import numpy as np
import matplotlib.pyplot as plt
from pandas.tools.plotting import autocorrelation_plot
from scipy import stats
import statsmodels.api as sm
import itertools

import biasWindow
biasWindow.windowSize = 10

extreme = 2147483648


def normalize(dy):

    while dy <= -extreme:
        dy += 2*extreme
        #print "Lower lim breached"
    while dy > extreme:
        dy -= 2*extreme
        #print "Higher lim breached"

    return dy


def getFiltered(g, filter):

    if type(g) == list:
        g = np.array(g)

    w = []
    a = []

    y = [g[:, 2], g[:, 6]]
    t = g[:, 0]/1.0e6

    dt = np.diff(t)

    clamp = np.vectorize(normalize)

    plt.close("all")

    for i in xrange(len(y)):
        w.append(np.zeros_like(t))
        a.append(np.zeros_like(t))

        dy = np.zeros_like(t)

        dy[1:] = clamp(np.diff(y[i]))/1.0e6

        w[i][1:] = np.divide(dy[1:], dt)

        if filter:

            w_raw = np.copy(w[i])
            w[i] =\
                   np.array(list(itertools.islice(biasWindow.centeredMean(w_raw),
                                                  w_raw.shape[0])))

            noise = w_raw - w[i]
            print "Noise for gyro ", (i+1)
            print "mu=", np.mean(noise), ",sig=", np.std(noise)

            plt.figure(i+1)
            plt.hist(noise, bins=200)
            plt.title("Noise distribution for gyro %d" % (i+1))

        a[i][1:] = np.divide(np.diff(w[i]), dt)

    plt.show()

    return t, y, w, a


def plot(g, PLOT):

    t, y, w, a = getFiltered(g, True)

    plt.close("all")

    for i in xrange(len(y)):

        print "Gryo # ", (i+1)

        print "Theta: mu=", np.min(y[i]), ",sig=", np.max(y[i])

        print "Omega: mu=", np.mean(w[i]), ",sig=", np.std(w[i])

        print "Alpha: mu=", np.mean(a[i]), ",sig=", np.std(a[i])

        if PLOT == 1:

            if False:
                plt.figure(i*3 + 1)
                plt.plot(t, y[i])
                plt.title("Yaw %d" % (i+1))
                plt.xlabel("Time (secs)")
                plt.ylabel("Yaw (degrees)")

            plt.figure(i*3 + 3)
            plt.plot(t, a[i])
            plt.title("Angular acceleration %d" % (i+1))
            plt.xlabel("Time (secs)")
            plt.ylabel("Angular acceleration (degrees/(sec^2))")

            plt.figure(i*3 + 2)
            plt.plot(t, w[i])
            plt.title("Yaw rate %d" % (i+1))
            plt.xlabel("Time (secs)")
            plt.ylabel("Yaw rate (degrees/sec)")

        elif PLOT == 2:

            plt.figure(i*1 + 1)
            plt.hist2d(w[i][:-1], a[i][1:], bins=150, normed=True)
            plt.title("Yaw rate change vs yaw rate for gyro %d" % (i+1))
            plt.xlabel("Yaw rate")
            plt.ylabel("Rate of change of yaw rate")
            plt.colorbar()

        elif PLOT == 3:

            plt.figure(i*1 + 1)
            plt.scatter(a[i][1:], w[i][:-1])
            plt.title("Yaw rate change vs yaw rate %d" % (i+1))
            plt.xlabel("Yaw rate")
            plt.ylabel("Rate of change of yaw rate")
            plt.grid(b=True)

    plt.show()


def regress(g, PLOT):

    t, y, w, a = getFiltered(g, True)

    for i in xrange(len(y)):

        print "For gyro ", (i+1)

        n = w[i].shape[0]
        mud = np.zeros(n)
        maxd = np.zeros(n)
        mind = np.zeros(n)
        muw = np.zeros(n)
        maxw = np.zeros(n)
        minw = np.zeros(n)
        d = np.zeros(n)

        for j in xrange(2, n):
            omega = w[i][max(0, j-5):j]
            alpha = a[i][max(0, j-5):j]

            mud[j] = np.mean(alpha)
            maxd[j] = np.max(alpha)
            mind[j] = np.min(alpha)

            muw[j] = np.mean(omega)
            maxw[j] = np.max(omega)
            minw[j] = np.min(omega)

            d[j] = omega[-1]-omega[0]

        model = sm.OLS(a[i][2:], np.column_stack((w[i][1:-1], mud[2:],
                                                  np.ones(w[i][1:-1].shape[0]))))
        results = model.fit()
        print results.summary()

        #print "Params are:", results.params
        print "Noise params (mu, sigma):", np.mean(results.resid), np.std(results.resid)

        if PLOT != 0:
            plt.close()

            plt.figure()
            plt.hist(results.resid, bins=200)
            plt.title("Residual distribution")

            plt.show()


def parse(files):

    PLOT = 1

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

                if count % 3.0e6 == 0:

                    print "Until observation ", count

                    regress(obs[j], 0)
                    #plot(obs[j], PLOT)

                    #return obs[j]

                    obs[j] = []

            #print "Until observation ", count
            #regress(obs[j])
            #plot(obs[j], PLOT)

            obs[j] = []

            print "\n"

    return obs

if __name__ == "__main__":
    files = ['../data/pimu_1.csv', '../data/pimu_3.csv']

    #files = ['../data/pimu_2.csv']

    obs = parse(files)
