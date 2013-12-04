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

    lim = None

    steps = None

    def __init__(self, obs, resW=0.0001, resB=0.001,
                 sigmaW=0.0085, sigmaB=0.00015, steps=1):
        """
        """
        print "FooBar"
        self.steps = steps

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
        lim = round(lim*1000.0)
        lim += (5-lim % 5)
        lim = lim/1000.

        self.states = np.concatenate((np.arange(0, -lim-self.resW,
                                                -self.resW)[:0:-1],
                                      np.arange(0, lim+self.resW, self.resW)),
                                     1)

        self.Nw = self.states.shape[0]

        self.Tw = np.empty(self.Nw)

        self.V = np.ones([self.N, self.Nw])*10.
        self.B = np.ones([self.N, self.Nw], dtype=np.int32)*-1

        self.Bp = np.zeros(self.N)

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

    def iterate(self, res, div, t, ol_idx, oh_idx, nl_idx, nh_idx, pml, mlidx):

        p_ml = pml
        ml = mlidx
        # looping over possible new states
        for i in xrange(nl_idx, nh_idx, (int)(res/self.resW)):
            p_max = -1e1000
            s_max = None
            # looping over old states
            for j in xrange(ol_idx, oh_idx, (int)(res/self.resW)):
                if self.Tw[abs(j-i)] <= 0 and self.V[t-self.steps, j] <= 0.:
                    p = self.V[t-self.steps, j] + self.Tw[abs(j-i)]
                    for sens in xrange(self.Ns):
                        oi = round(self.obs[sens, t-self.steps]*div)
                        of = round(self.obs[sens, t]*div)

                        if res == 0.005:
                            r = oi % 5
                            oi += 5-r if r > 2 else -r

                            r = of % 5
                            of += 5-r if r > 2 else -r

                            bi = (oi/div) - round(self.states[j], 3)
                            bf = (of/div) - round(self.states[i], 3)

                        else:
                            bi = (oi/div) - round(self.states[j], 4)
                            bf = (of/div) - round(self.states[i], 4)

                        p += self.getBiasTrans(bi, bf)

                    if p > p_max:
                        p_max = p
                        s_max = j

            self.V[t, i] = p_max
            self.B[t, i] = s_max

            if p_max > p_ml:
                p_ml = p_max
                ml = i

        return ml

    def findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """
        for i in xrange(self.Nw):
            p = self.getProb(self.states[i], 0, 0.5, self.resW)
            if p != 0:
                self.V[0, i] = math.log(p)
            else:
                self.V[0, i] = -1e1000

        self.B[0] = np.arange(self.Nw, dtype=np.int32)
        self.Bp[0] = 0.

        for t in xrange(self.steps, self.N, self.steps):  # looping over time
            div = 1000.0
            ml = self.iterate(0.005, div, t, 0, self.Nw,
                              0, self.Nw, -1e1000, None)

            self.Bp[t] = self.states[ml]

            refine = True

            if refine:
                src = self.B[t, ml]
                ol_idx = max(src, 0)
                oh_idx = min(src+1, self.Nw)

                if self.V[t, ol_idx] == 10.:
                    ol_idx = src
                if self.V[t, oh_idx-1] == 10.:
                    oh_idx = src+1

                nl_idx = max(ml - (int)(0.005/self.resW)+1, 0)
                nh_idx = min(ml + (int)(0.005/self.resW), self.Nw)

                div = 10000.
                ml_fine = self.iterate(0.0001, div, t, ol_idx, oh_idx,
                                       nl_idx, nh_idx, self.V[t, ml], ml)

                self.Bp[t] = self.states[ml_fine]

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
