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

    def __init__(self, obs, resW=0.005, resB=0.001,
                 sigmaW=0.0085, sigmaB=0.00025):
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

        self.states = np.concatenate((np.arange(0, -lim, -self.resW)[:0:-1],
                                      np.arange(0, lim, self.resW)), 1)

        self.Nw = self.states.shape[0]

        self.Tw = np.empty(self.Nw)

        self.V = np.ones([self.N, self.Nw])*10
        self.B = np.empty([self.N, self.Nw])

        self.Bp = np.empty(self.N)

    def createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """

        prob = stats.norm(0, self.sigmaW)

        for i in xrange(self.Nw):
            p = prob.pdf(i*self.resW)/100.0
            if p != 0:
                self.Tw[i] = math.log(p)
            else:
                self.Tw[i] = 10.

    def getProb(self, x, mu, sigma, res):
        return (0.5*(1+math.erf((x-mu+res)/(math.sqrt(2)*sigma))) -
                0.5*(1+math.erf((x-mu-res)/(math.sqrt(2)*sigma))))

    def getBiasTrans(self, init, final):
        """
        """
        delB = final - init
        prob = self.getProb(delB, 0, self.sigmaB, 0.00001)
        if prob != 0:
            return math.log(prob)
        else:
            return -9999

    def findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """
        for i in xrange(self.Nw):
            self.V[0, i] = math.log(self.getProb(self.states[i], 0,
                                                 0.5, self.resW))

        self.B[0] = np.arange(self.Nw)
        self.Bp[0] = 0.

        div = 1000.0

        for t in xrange(1, self.N):  # looping over time
            p_ml = -1e1000
            ml = None
            for i in xrange(self.Nw):  # looping over possible new states
                p_max = -1e1000
                s_max = None
                # looping over old states
                for j in xrange(self.Nw):
                    if self.Tw[abs(j-i)] <= 0.:
                        p = self.V[t-1, j] + self.Tw[abs(j-i)]
                        for sens in xrange(self.Ns):
                            oi = round(self.obs[sens, t-1]*div)
                            of = round(self.obs[sens, t]*div)

                            r = oi % 5
                            oi += 5-r if r > 2 else -r

                            r = of % 5
                            of += 5-r if r > 2 else -r

                            bi = (oi/div) - round(self.states[j], 3)
                            bf = (of/div) - round(self.states[i], 3)
                            p += self.getBiasTrans(bi, bf)

                        if p > p_max:
                            p_max = p
                            s_max = j

                self.V[t, i] = p_max
                self.B[t, i] = s_max

                if p_max > p_ml:
                    p_ml = p_max
                    ml = i
            self.Bp[t] = self.states[ml]

        #st = np.argmax(self.V[-1])

        #for t in xrange(self.N-1, -1, -1):
        #    self.Bp[t] = self.states[st]
        #    st = self.B[t, st]

    def run(self):
        """
        Run the viterbi algorithm, and plot data against GT
        """
        self.createMatrices()
        self.findSequence()

        plt.plot(self.Bp, 'r')
        plt.plot(self.omega, 'g')
# if __name__ == "__main__":
# bias, omega, obs = generateData.generate(2, 150000)
