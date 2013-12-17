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
    Bpr = None

    lim = None

    Vs = None

    def __init__(self, obs, resW=0.005, resB=0.001,
                 sigmaW=0.0085, sigmaB=0.0015):
        """
        """
        self.obs = obs

        self.resB = resB
        self.resW = resW

        self.sigmaW = sigmaW
        self.sigmaB = sigmaB

        self.Tb = stats.norm(0, sigmaB)

        self.N = self.obs.shape[1]
        self.Ns = self.obs.shape[0]

        minw = 1e10
        maxw = -1e10

        for i in xrange(self.Ns):
            val = np.amin(self.obs[i, :])
            if val < minw:
                minw = val

            val = np.amax(self.obs[i, :])
            if val > maxw:
                maxw = val

        lim = max((minw*-1), maxw)
        self.lim = lim

        self.states = np.concatenate((np.arange(0, -lim,
                                                -self.resW)[:0:-1], np.arange(0, lim, self.resW)), 1)

        self.Nw = self.states.shape[0]

        self.Tw = np.ones([self.Nw, self.Nw])*10.

        self.V = np.ones([self.N, self.Nw])*10.
        self.Vs = np.ones([2, self.Nw])*10.
        self.B = np.ones([self.N, self.Nw], dtype=np.int32)*-1

        self.Bp = np.zeros(self.N)
        self.Bpr = np.zeros(self.N)

    def createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """

        prob = stats.norm(0, self.sigmaW)

        for j in xrange(self.Nw):
            for i in xrange(self.Nw):
                # probability that state j came from state i
                self.Tw[i, j] = prob.pdf(abs(j-i)*self.resW)/100.

            #self.Tw[:, j] /= sum(self.Tw[:, j])

            for i in xrange(self.Nw):
                p = self.Tw[i, j]
                if p != 0:
                    self.Tw[i, j] = math.log(p)
                else:
                    self.Tw[i, j] = 10.

    def getProb(self, x, mu, sigma, res):
        return (0.5*(1+math.erf((x-mu+res)/(math.sqrt(2)*sigma))) -
                0.5*(1+math.erf((x-mu-res)/(math.sqrt(2)*sigma))))

    def getPDF(self, x, mu, sigma):
        return (1/(math.sqrt(2)*sigma))*math.exp((-0.5)*pow(((x-mu)/sigma), 2))

    def getBiasTrans(self, init, final):
        """
        """
        delB = final - init
        prob = self.getProb(delB, 0, self.sigmaB, 0.00001)
        if prob != 0:
            return math.log(prob)
        else:
            return -9999

    def iterate(self, res, t):

        p_ml = -1e10
        ml = None

        # looping over possible new states
        for i in xrange(self.Nw):
            p_max = -1e10
            s_max = None
            # looping over old states
            for j in xrange(self.Nw):
                if self.Tw[j, i] <= 0. and self.Vs[0, j] <= 0.:
                    p = self.Vs[0, j] + self.Tw[j, i]
                    for sens in xrange(self.Ns):
                        bi = self.obs[sens, t-1] - self.states[j]
                        bf = self.obs[sens, t] - self.states[i]

                        p += self.getBiasTrans(bi, bf)

                    if p > p_max:
                        p_max = p
                        s_max = j

            self.Vs[1, i] = p_max
            self.B[t, i] = s_max

            if p_max > p_ml:
                p_ml = p_max
                ml = i
                self.Bp[t] = self.states[ml]

        self.Vs[0] = np.copy(self.Vs[1])

        return ml

    def findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """
        for i in xrange(self.Nw):
            p = self.getProb(self.states[i], 0., 0.5, self.resW)
            if p != 0:
                self.Vs[0, i] = math.log(p)
            else:
                self.Vs[0, i] = 1e2

        self.B[0] = np.arange(self.Nw, dtype=np.int32)
        self.Bp[0] = 0.

        for t in xrange(1, self.N):  # looping over time
            ml = self.iterate(0.005, t)

        st = ml

        for t in xrange(self.N-1, -1, -1):
            self.Bpr[t] = self.states[st]
            st = self.B[t, st]

    def createRun(self, omega):
        """
        Run the viterbi algorithm, and plot data against GT
        """
        self.createMatrices()
        self.findSequence()

        plt.plot(self.Bpr, 'r')
        plt.plot(omega[:self.Bp.shape[0]], 'g')
        plt.show()

    def run(self, omega):

        self.findSequence()

        plt.plot(self.Bp, 'r')
        plt.plot(omega[:self.Bp.shape[0]], 'g')

        plt.show()
