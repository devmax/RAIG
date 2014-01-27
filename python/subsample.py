import csv
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import itertools

extreme = 2147483648
window = 11


def normalize(dy):

    while dy <= -extreme:
        dy += 2*extreme
        #print "Lower lim breached"
    while dy > extreme:
        dy -= 2*extreme
        #print "Higher lim breached"

    return dy


# subsample the omega readings over a window size, and return the new
# subsampled omega as well as its diff array

def getFiltered(g, PLOT):

    if type(g) == list:
        g = np.array(g)

    w = []
    a = []

    y = [g[:, 2], g[:, 6]]
    t = g[:, 0]/1.0e6

    subsize = int(t.shape[0]/window)

    t = t[:subsize*window]

    dt = np.diff(t)

    clamp = np.vectorize(normalize)

    plt.close("all")

    for i in xrange(len(y)):

        y[i] = y[i][:subsize*window]

        w.append(np.zeros(subsize*window))
        a.append(np.zeros(subsize))

        dy = np.zeros_like(y[i])

        dy[1:] = clamp(np.diff(y[i]))/1.0e6

        w[i][1:] = np.divide(dy[1:], dt)

        w_raw = w[i][window/2::window]
        w[i] = np.mean(w[i].reshape(-1, window), 1)
        y[i] = np.mean(y[i].reshape(-1, window), 1)

        a[i][1:] = np.diff(w[i])  # acceleration here is not divided
                                  # by dt, to make the regression well
                                  # conditioned

        noise = w_raw - w[i]

        if PLOT != 0:

            plt.figure()
            plt.hist(noise, bins=250)
            plt.title("Gaussian Noise distribution")

        print "Gaussian noise:(mu, sig)=", np.mean(noise), np.std(noise)

    plt.show()

    t = np.mean(t.reshape(-1, window), 1)

    return t, y, w, a


def plot(t, y, w, a, PLOT):

    plt.close("all")

    for i in xrange(len(w)):

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


def regress(subsampled, PLOT):

    t, w, a = subsampled[0], subsampled[1], subsampled[2]

    for i in xrange(len(w)):

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

        model = sm.OLS(a[i][2:], np.column_stack((w[i][1:-1],# mud[2:],
                                                  np.ones(w[i][1:-1].shape[0]))))
        results = model.fit()
        print results.summary()

        #print "Params are:", results.params
        print "Residual noise params (mu, sigma):", np.mean(results.resid), np.std(results.resid)

        if PLOT != 0:
            plt.close()

            plt.figure()
            plt.hist(results.resid, bins=200)
            plt.title("Residual distribution")

            plt.show()


def parse(files):

    PLOT = 0

    obs = []

    subsampled = [list() for i in xrange(len(files))]

    for j in xrange(len(files)):

        print "Gyroscope ", j+1
        print "*"*50

        with open(files[j]) as datafile:

            subsampled[j].append(np.array([]))  # t
            subsampled[j].append([np.array([]), np.array([])])  # w
            subsampled[j].append([np.array([]), np.array([])])  # a

            data = csv.reader(datafile)
            data.next()

            count = 0

            for i in xrange(int(20.e3)):
                data.next()

            for row in data:
                obs.append([float(item) for item in row])

                count += 1
                if count > 1.0e7:
                    break

                if count % int(window*2.0e5) == 0 or count >= 1.0e7:

                    print "Until observation ", count

                    t, y, w, a = getFiltered(obs, PLOT)

                    subsampled[j][0] = np.append(subsampled[j][0], t)
                    subsampled[j][1][0] = np.append(subsampled[j][1][0], w[0])
                    subsampled[j][1][1] = np.append(subsampled[j][1][1], w[1])
                    subsampled[j][2][0] = np.append(subsampled[j][2][0], a[0])
                    subsampled[j][2][1] = np.append(subsampled[j][2][1], a[1])

                    obs = []

            #plot(subsampled, PLOT)
            regress(subsampled[j], PLOT)

            print "\n"

    return subsampled

if __name__ == "__main__":
    files = ['../data/pimu_1.csv', '../data/pimu_3.csv']

    #files = ['../data/pimu_2.csv']

    obs = parse(files)
