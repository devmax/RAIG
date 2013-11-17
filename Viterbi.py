import numpy as np
import scipy.stats as stats
import math

import matplotlib.pyplot as plt
# import generateData


class Viterbi():
    """
    """
    obs = None
    resW = None
    resB = None

    sigmaW = None
    sigmaB = None

    N = None
    Nw = None
    Nb = None
    Ns = None

    states = None
    biases = None

    V = None
    B = None

    Bp = None

    def __init__(self, obs, resW=0.05, resB=0.001,
                 sigmaW=0.0085, sigmaB=0.00025):
        """
        """
        self.obs = obs
        self.resW = resW
        self.resB = resB

        self.sigmaW = sigmaW
        self.sigmaB = sigmaB

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

        self.states = np.arange(-lim, lim, self.resW)
        self.biases = np.arange(-0.5, 0.5, self.resB)

        self.Nw = self.states.shape[0]
        self.Nb = self.biases.shape[0]

        self.Tw = np.empty([self.Nw, self.Nw])
        self.Tb = np.empty([self.Nb, self.Nb])

        self.V = np.empty([self.Ns, self.N])
        self.B = np.empty([self.Ns, self.N])

        self.Bp = np.empty([1, self.N])

    def createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """
        prob = stats.norm(0, self.sigmaW)

        for i in xrange(self.Nw):
            for j in xrange(self.Nw):
                self.Tw[i, j] = prob.pdf((j-i)*self.resW)
            self.Tw[i, :] /= np.sum(self.Tw[i, :])

        for i in xrange(self.Nb):
            for j in xrange(self.Nb):
                self.Tb[i, j] = prob.pdf((j-i)*self.resB)
            self.Tb[i, :] /= np.sum(self.Tb[i, :])

    def getBiasTrans(self, init, final):
        """
        """
        fr = (init - self.biases[0])/self.resB
        to = (final - self.biases[0])/self.resB

        return self.Tb[fr, to]

    def findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """

        bias = np.empty([self.Ns, self.N])

        priorState = stats.norm(0, 0.75)

        self.V[0, 0] = math.log(priorState.cdf(self.states[0]))
        self.V[0, self.Nw-1] = math.log(1 - priorState.cdf(self.states[-1]))

        self.B[0, 0] = 0
        self.B[0, self.Nw-1] = self.Nw-1

        for i in xrange(1, self.Nw-1):
            self.V[0, i] = math.log(priorState.cdf(self.states[i+1]) -
                                    priorState.cdf(self.states[i-1]))
            for j in xrange(self.Ns):
                bias[j, 0] = self.obs[j, 0]
                self.B[0, i] = i

        for t in xrange(1, self.N):  # looping over time
            for i in xrange(self.Nw):  # looping over possible new states
                p_max = -1
                s_max = None
                # looping over old states
                for j in xrange(self.Nw):
                    p = self.V[t-1, j] + math.log(self.Tw[j, i])
                    for sens in xrange(self.Ns):
                        bi = self.obs[sens, t-1] - self.states[j]
                        bf = self.obs[sens, t] - self.states[i]
                        if bf < self.biases[0] or bf > self.biases[-1]:
                            raise Exception("Bias has exceeded limit!")
                        p += math.log(self.getBiasTrans(bi, bf))

                    if p > p_max:
                        p_max = p
                        s_max = j

                self.V[t, i] = p_max
                self.B[t, i] = s_max

        st = np.argmax(self.V[-1, :])

        for t in xrange(self.N-1, -1, -1):
            self.Bp[t] = self.states[st]
            st = self.B[t, st]

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
