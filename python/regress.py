import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

import parse


def regress(data, PLOT):

    for i in xrange(len(data)):

        dt, w, dw = data[i][:, 0], data[i][:, 1], data[i][:, 2]

        print "For gyro ", (i+1)

        model = sm.OLS(dw, np.ones(w.shape[0]))  # np.column_stack((w, np.sign(w))))
        results = model.fit()

        print results.summary()

        #print "Params are:", results.params
        print "Noise params (mu, sigma):", np.mean(results.resid), np.std(results.resid)

        if PLOT != 0:
            plt.close("all")

            plt.figure()
            plt.hist(results.resid, bins=200)
            plt.title("Residual distribution")

            plt.show()

if __name__ == "__main__":

    files = ['/home/dev/Documents/RAIG/data/test1_2_reg.txt',
             '/home/dev/Documents/RAIG/data/test2_2_reg.txt',
             '/home/dev/Documents/RAIG/data/test3_2_reg.txt']

    data = parse.separate(files)

    regress(data, 0)
