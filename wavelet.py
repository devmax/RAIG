import pywt
import matplotlib.pyplot as plt


def getFeatureError(data, gtruth, start, end, family):

    obs = data[:, start:end]
    omega = gtruth[start:end]

    num = obs.shape[0]
    coeffs = []
    thresh = 1e-4

    for i in xrange(num):
        coeffs.append(pywt.wavedec(obs[i], family))
        #print "For sensor ", i+1
        #print coeffs[i]

    for lev in xrange(len(coeffs[0])):
        for j in xrange(len(coeffs[0][lev])):
            mean = 0
            err = 0
            for k in xrange(num):
                mean += coeffs[k][lev][j]
            mean /= num

            for k in xrange(num):
                err += (coeffs[k][lev][j] - mean)**2
            if mean != 0:
                err = abs(err/mean)
            else:
                err = 0

            if (err > thresh):
                for k in xrange(num):
                    coeffs[k][lev][j] = 0

    plt.figure()
    plt.hold(True)

    plt.plot(omega, 'g')

    for i in xrange(num):
        rec = pywt.waverec(coeffs[i], family)
        plt.plot(rec)

    plt.show()
