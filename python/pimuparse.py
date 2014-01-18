import csv
import numpy as np
import matplotlib.pyplot as plt
from pandas.tools.plotting import autocorrelation_plot
from scipy import stats
import statsmodels.api as sm

import biasWindow

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

    y1 = (g[:, 2])
    y2 = (g[:, 6])
    t = g[:, 0]/1.0e6

    w1 = np.zeros_like(t)
    w2 = np.zeros_like(t)

    a1 = np.zeros_like(t)
    a2 = np.zeros_like(t)

    for i in xrange(1, t.shape[0]):
        dt = (t[i]-t[i-1])
        dy = [y1[i]-y1[i-1], y2[i]-y2[i-1]]

        for j in xrange(len(dy)):
            dy[j] = normalize(dy[j])/1.0e6

        w1[i] = dy[0]/dt
        w2[i] = dy[1]/dt

    if filter:
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
    else:
        for i in xrange(1, t.shape[0]):
            dt = (t[i]-t[i-1])
            a1[i] = (w1[i]-w1[i-1])/dt
            a2[i] = (w2[i]-w2[i-1])/dt

    return t, y1, y2, w1, w2, a1, a2


def compareLL(g):

    t, y1, y2, w1, w2, a1, a2 = getFiltered(g, False)

    t1 = 50.0
    t2 = 20.0
    diff1 = [0, 0]
    diff2 = [0, 0]

    for i in xrange(w1.shape[0]):
        diff1[0] += (w1[i]*t1 - a1[i])**2
        diff1[1] += (w2[i]*t1 - a2[i])**2

        diff2[0] += (w1[i]*t2 - a1[i])**2
        diff2[1] += (w2[i]*t2 - a2[i])**2

    print "SSE for theta=50: [", diff1[0], ",", diff1[1], "]"
    print "SSE for theta=20: [", diff2[0], ",", diff2[1], "]"


def plot(g, PLOT):

    t, y1, y2, w1, w2, a1, a2 = getFiltered(g, False)

    print "Theta: [", np.min(y1), ",", np.max(y1), "]", "[",
    print np.min(y2), ",", np.max(y2), "]"

    print "Omega: [", np.mean(w1), ",", np.std(w1), "]", " [",
    print np.mean(w2), ",", np.std(w2), "]\n"

    print "Alpha: [", np.mean(a1), ",", np.std(a1), "]", " [",
    print np.mean(a2), ",", np.std(a2), "]\n"

    plt.close()

    if PLOT == 1:

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

    elif PLOT == 2:

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

    elif PLOT == 3:

        plt.close()

        plt.figure(1)
        plt.scatter(w1[:-1], w1[1:])
        plt.title("Yaw rate change vs yaw rate(1)")
        plt.xlabel("Yaw rate")
        plt.ylabel("Rate of change of yaw rate")
        plt.grid(b=True)

        plt.figure(2)
        plt.scatter(w2[:-1], w2[1:])
        plt.title("Yaw rate change vs yaw rate(2)")
        plt.xlabel("Yaw rate")
        plt.ylabel("Rate of change of yaw rate")
        plt.grid(b=True)

    plt.show()


def regress(g):

    t, y1, y2, w1, w2, a1, a2 = getFiltered(g, False)

    print "First axis:"

    model = sm.OLS(w1[1:], np.column_stack((w1[:-1], np.ones(a1.shape[0]-1))))
    results = model.fit()

    print "Params are:", results.params
    print "Noise params (mu, sigma):", np.mean(results.resid), np.std(results.resid)
    plt.figure()
    plt.hist(results.resid, bins=200)
    plt.title("Residual distribution")

    #print results.summary()

    print "Second axis:"

    model = sm.OLS(w2[1:], np.column_stack((w2[:-1], np.ones(a2.shape[0]-1))))
    results = model.fit()

    print "Params are:", results.params
    print "Noise params (mu, sigma):", np.mean(results.resid), np.std(results.resid)
    plt.figure()
    plt.hist(results.resid, bins=200)
    plt.title("Residual distribution")
    #print results.summary()

    plt.show()


def parse(files):

    PLOT = 3

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

                if count % 3.e6 == 0:

                    print "Until observation ", count

                    regress(obs[j])
                    #plot(obs[j], PLOT)
                    #compareLL(obs[j])

                    #return obs[j]

                    obs[j] = []

            #print "Until observation ", count
            #regress(obs[j])
            #plot(obs[j], PLOT)
            #compareLL(obs[j])

            obs[j] = []

            print "\n"

    return obs

if __name__ == "__main__":
    files = ['../data/pimu_2.csv', '../data/pimu_1.csv',
             '../data/pimu_3.csv']

    #files = ['../data/pimu_2.csv']

    obs = parse(files)
