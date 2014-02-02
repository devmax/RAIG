import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import math


def logProb(x, y, w, sigma):

    return -(0.5*math.log(2*np.pi*(sigma**2)) + 0.5*((y-np.dot(w, x))/sigma)**2)


if __name__ == "__main__":

    x = np.arange(500)
    y = np.arange(150, 650) #+ np.random.normal(0, 1.5, 500)

    x = np.column_stack((x, np.ones(y.shape[0])))

    model = sm.OLS(y, x)
    res = model.fit()

    ll = 0.0
    sigma = np.std(res.resid)

    for i in xrange(y.shape[0]):
        ll += logProb(x[i], y[i], res.params, sigma)

    print "My LL: ", ll

    print res.summary()
