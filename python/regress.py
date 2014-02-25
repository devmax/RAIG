import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

import parse


def regress(data, PLOT):

    for i in xrange(len(data)):

        dt, w, dw = data[i][:, 0], data[i][:, 1], data[i][:, 2]

        print "For gyro ", (i+1)

<<<<<<< HEAD
        model = sm.OLS(dw, np.ones(w.shape[0]))  # np.column_stack((w, np.sign(w))))
=======
        model = sm.OLS(dw, np.column_stack((w, np.sign(w))))  # np.ones(w.shape[0]))))
>>>>>>> 4664a91daa215b619bebb03d44d5a1e894575ba0
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

<<<<<<< HEAD
    files = ['/home/dev/Documents/RAIG/data/test1_2_reg.txt',
             '/home/dev/Documents/RAIG/data/test2_2_reg.txt',
             '/home/dev/Documents/RAIG/data/test3_2_reg.txt']
=======
    files = ['../data/big0_reg.txt', '../data/big1_reg.txt',
             '../data/big2_reg.txt']
>>>>>>> 4664a91daa215b619bebb03d44d5a1e894575ba0

    data = parse.separate(files)

    regress(data, 0)
