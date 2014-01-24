import numpy as np
import scipy.stats as stats
import math

import matplotlib.pyplot as plt


class Viterbi():
    """
    """
    obs = None
    resW = None
    resB = None

    sigmaW = None
    sigmaB = None

    Tw = None
    Tb = None

    N = None
    Nw = None
    Ns = None

    states = None

    V = None
    B = None

    Bp = None
    Bpv = None

    lim = None

    Vs = None

    def __init__(self, obs, omega, resW=0.005, resB=0.001,
                 sigmaW=0.0085, sigmaB=0.00015):
        """
        """
        self.obs = obs

        self.omega = omega[:obs.shape[1]]

        self.resB = resB
        self.resW = resW

        self.sigmaW = sigmaW
        self.sigmaB = sigmaB

        self.theta = [0., -0.01, 0.75]

        self.Tb = stats.norm(0, sigmaB)

        self.N = self.obs.shape[1]
        self.Ns = self.obs.shape[0]

        minw = np.amin(self.obs)
        maxw = np.amax(self.obs)

        lim = max((minw*-1), maxw)
        self.lim = lim

        self.states = np.concatenate((np.arange(0, -lim,
                                                -self.resW)[:0:-1], np.arange(0, lim, self.resW)), 1)

        self.Nw = self.states.shape[0]

        self.Vs = np.zeros([2, self.Nw])

        self.Tw = np.zeros([self.Nw, self.Nw])

        self.reinitialize()
        self.createMatrices()

    def createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """
        for j in xrange(self.Nw):
            lp = np.vectorize(self.logProb)

            self.Tw[:, j] = lp(np.abs(self.states[j]-self.states), 0,
                               self.sigmaW)

    def getProb(self, x, mu, sigma, res):
        return (0.5*(1+math.erf((x-mu+res)/(math.sqrt(2)*sigma))) -
                0.5*(1+math.erf((x-mu-res)/(math.sqrt(2)*sigma))))

    def logProb(self, x, mu, sigma):

        return -(np.log(sigma)+0.5*((x-mu)/sigma)**2)

    def getBiasTrans(self, init, final, mu):
        """
        """
        #delB = final - init
        #prob = self.getProb(delB, 0, self.sigmaB, 0.00001)

        return self.logProb(final-init, self.theta[0] +
                            self.theta[1]*init + self.theta[2]*mu, 0.0068)

    def iterate(self, t):
        p_ml = -1e10
        ml = None

        bias = np.vectorize(self.getBiasTrans)

        # looping over possible new states
        for i in xrange(self.Nw):
            p_max = -1e10
            s_max = None

            p = self.Vs[0] + self.Tw[:, i]

            for sens in xrange(self.Ns):
                bi = self.obs[sens, t-1]-self.states
                bf = self.obs[sens, t]-self.states[i]

                p += bias(bi, bf)

            p_max = np.max(p)
            s_max = np.argmax(p)

            if s_max is None:
                print "(t, i) is ", t, i

            self.Vs[1, i] = p_max
            self.B[t, i] = s_max

            if p_max > p_ml:
                p_ml = p_max
                ml = i
                self.Bp[t] = self.states[ml]

        self.Vs[0] = np.copy(self.Vs[1])

        return ml

    def reinitialize(self):

        lp = np.vectorize(self.logProb)
        self.Vs[0] = lp(self.states, 0., 0.5)

        self.Vs[1] = np.ones(self.Vs.shape[1])*100.
        self.B = np.ones([self.N, self.Nw], dtype=np.int32)*-1
        self.B[0] = np.arange(self.Nw, dtype=np.int32)
        self.Bp = np.zeros(self.N)
        self.Bpv = np.zeros(self.N)

    def findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """
        self.reinitialize()

        for t in xrange(1, self.N):  # looping over time
            ml = self.iterate(t)

        st = ml

        for t in xrange(self.N-1, -1, -1):
            self.Bpv[t] = self.states[st]
            st = self.B[t, st]

    def plot(self):
        """
        Run the viterbi algorithm, and plot data against GT
        """

        plt.subplot(211)

        plt.plot(self.Bp, 'r', label="Estimate")
        plt.plot(self.omega[:self.Bp.shape[0]], 'g', label="Ground truth")
        plt.xlabel("Time")
        plt.ylabel("Rate")

        plt.subplot(212)
        plt.plot(self.Bp-self.omega[:self.Bp.shape[0]], 'b', label="Error")
        plt.xlabel("Time")
        plt.ylabel("Rate Error")

        plt.show()

    def run(self):

        self.findSequence()

        self.plot()
