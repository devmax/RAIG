import numpy as np
import parse

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

    def __init__(self, omega, N=100, Ns=2, resW=1, sigmaW=0.1,
                 sigmaB=7.0):
        """
        """

        self.scale = 131.0

        self.N = N
        self.Ns = Ns

        self.b = parse.separate(maxCount=self.N, maxFiles=self.Ns)

        assert(self.Ns == len(self.b))

        self.omega = omega[:self.N]

        self.obs = np.zeros([self.Ns, self.N])

        for i in xrange(self.Ns):
            self.obs[i] = self.b[i][:, 1] + self.omega

        self.resW = resW

        self.sigmaW = sigmaW
        self.sigmaB = sigmaB

        self.theta = [0.]

        self.Ns = self.obs.shape[0]

        minw = np.amin(self.obs)
        maxw = np.amax(self.obs)

        lim = max((minw*-1), maxw)
        self.lim = lim

        self.states = np.concatenate((np.arange(0, -lim,
                                                -self.resW)[:0:-1], np.arange(0, lim, self.resW)), 1)

        self.Nw = self.states.shape[0]

        self.Vs = np.zeros([2, self.Nw])

        self.Tw = np.zeros(self.Nw)

        self.reinitialize()
        self.createMatrices()

    def createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """

        for i in xrange(self.Nw):
            self.Tw[i] = self.logProb(i*self.resW, 0.0, round(self.sigmaW*self.scale))

        if False:
            for j in xrange(self.Nw):
                lp = np.vectorize(self.logProb)

                self.Tw[:, j] = lp(np.abs(self.states[j]-self.states),
                                   0, self.sigmaW)

    def logProb(self, x, mu, sigma):

        return -(np.log(sigma)+0.5*((x-mu)/sigma)**2)

    def getBiasTrans(self, init, final):
        """
        """
        return self.logProb(final-init, self.theta[0], self.sigmaB)

    def iterate(self, t):
        p_ml = -1e10
        ml = None

        allp = []

        # looping over possible new states
        for i in xrange(self.Nw):
            p_max = -1e10
            s_max = None

            plist = []

            # loop over old states that were possible sources
            for j in xrange(self.Nw):
                p = self.Vs[0, j] + self.Tw[abs(j-i)]

                for sens in xrange(self.Ns):
                    bi = self.obs[sens, t-1]-self.states[j]
                    bf = self.obs[sens, t]-self.states[i]

                    p += self.getBiasTrans(bi, bf)

                plist.append(p)

                if p > p_max:
                    p_max = p
                    s_max = j

            allp.append(np.array(plist))

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
        self.Vs[0] = lp(self.states, 0., 13.0)

        self.Vs[1] = np.ones(self.Vs.shape[1])
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
        plt.figure(1)
        plt.subplot(211)

        plt.plot(self.Bp, 'r', label="Estimate")
        plt.plot(self.omega[:self.Bp.shape[0]], 'g', label="Ground truth")
        plt.xlabel("Time")
        plt.ylabel("Rate")

        plt.legend()

        plt.subplot(212)
        plt.plot(self.Bp-self.omega[:self.Bp.shape[0]], 'b', label="Error")
        plt.xlabel("Time")
        plt.ylabel("Rate Error")

        plt.legend()

        plt.show()

    def run(self):

        self.findSequence()

        self.plot()
