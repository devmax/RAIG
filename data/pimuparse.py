import csv
import numpy as np
import matplotlib.pyplot as plt


def parse(files):

    obs = [list() for i in xrange(len(files))]

    for j in xrange(len(files)):
        with open(files[j]) as datafile:
            data = csv.reader(datafile)
            data.next()
            for i in xrange(int(10.5e6)):
                data.next()

            for row in data:
                obs[j].append([float(item) for item in row])

        g = np.array(obs[j])

        y1 = g[:, 2]/1.0e6
        y2 = g[:, 6]/1.0e6
        t = g[:, 0]/1.0e6
        t = t - t[0]

        dw1 = np.zeros_like(t)
        dw2 = np.zeros_like(t)

        da1 = np.zeros_like(t)
        da2 = np.zeros_like(t)

        for i in xrange(1, t.shape[0]):
            dt = (t[i]-t[i-1])
            dw1[i] = (y1[i]-y1[i-1])/dt
            dw2[i] = (y2[i]-y2[i-1])/dt

            da1[i] = (dw1[i]-dw1[i-1])/dt
            da2[i] = (dw2[i]-dw2[i-1])/dt

        plt.figure()
        plt.plot(t, y1)

        plt.figure()
        plt.plot(t, y2)

        plt.figure()
        plt.plot(t, dw1)

        plt.figure()
        plt.plot(t, dw2)

        plt.figure()
        plt.plot(t, da1)

        plt.figure()
        plt.plot(t, da2)

        plt.show()

    return obs

if __name__ == "__main__":
    files = ['pimu_1.csv', 'pimu_2.csv',
             'pimu_3.csv']

    obs = parse(files)
