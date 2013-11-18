import numpy as np
import scipy.stats as stats
import math
cimport numpy as np

import matplotlib.pyplot as plt

cdef class Viterbi:
    """
    """
    DTYPE = np.float64
    ctypedef np.float64_t DTYPE_t

    cdef double resW
    cdef double resB

    cdef double sigmaW
    cdef double sigmaB

    cdef double Tw
    cdef double Tb

    cdef unsigned int N
    cdef unsigned int Nw
    cdef unsigned int Nb
    cdef unsigned int Ns

    cdef np.ndarray[DTYPE_t, ndim=1] states

    cdef np.ndarray[DTYPE_t, ndim=2] obs
    cdef np.ndarray[DTYPE_t, ndim=2] V
    cdef np.ndarray[DTYPE_t, ndim=2] B

    cdef np.ndarray[DTYPE_t, ndim=1] Bp
    
    cdef double sqrt_tpi
    cdef double sqrt_two

    cdef inline double double_max(double a, double b): return a if a>=b else b
    cdef inline double double_min(double a, double b): return a if a<=b else b

    cdef extern from "math.h":
        double exp(double x)
        double log(double x)
        double erf(double x)
        double sqrt(double x)
        double pow(double base, double exponent)

    def __init__(self,np.ndarray[DTYPE_t, ndim=2] obs, double resW=0.005, double resB=0.001,
                 double sigmaW=0.0085, double sigmaB=0.00025):
        """
        """
        self.sqrt_tpi = sqrt(2*math.pi)
        self.sqrt_two = sqrt(2)

        self.obs = obs

        self.resB = resB
        self.resW = resW

        self.sigmaW = sigmaW
        self.sigmaB = sigmaB

        self.N = self.obs.shape[1]
        self.Ns = self.obs.shape[0]

        cdef double minw = 1e10
        cdef double maxw = -1e10

        cdef double val
        for i in xrange(self.Ns):
            val = np.amin(self.obs[i, :])
            if val < minw:
                minw = val

            val = np.amax(self.obs[i, :])
            if val > maxw:
                maxw = val

        cdef double lim = double_max((minw*-1), maxw)

        self.states = np.concatenate((np.arange(-lim, 0, self.resW),
                                      np.arange(0, lim, self.resW)), 1)

        self.Nw = self.states.shape[0]

        self.Tw = np.empty([self.Nw, self.Nw])

        self.V = np.empty([self.N, self.Nw])
        self.B = np.empty([self.N, self.Nw])

        self.V[:, :] = -15000

        self.Bp = np.empty(self.N)

    cpdef double getPDF(self, double x, double mu, double sigma):
        return (1/(self.sqrt_tpi*sigma))*exp((-0.5)*pow(((x-mu)/sigma), 2))

    cpdef double getProb(self, double x, double mu, double sigma, double res):
        return (0.5*(1+erf((x+res-mu)/(self.sqrt_two*sigma))) -
                0.5*(1+erf((x-res-mu)/(self.sqrt_two*sigma))))

    cpdef createMatrices(self):
        """
        Create state transition matrices, and observation matrices
        """
        for i in range(self.Nw):
            for j in range(self.Nw):
                self.Tw[i, j] = getPDF((j-i)*self.resW, 0, self.sigmaW)
            self.Tw[i, :] /= np.sum(self.Tw[i, :])
            for j in range(self.Nw):
                if self.Tw[i, j] != 0:
                    self.Tw[i, j] = log(self.Tw[i, j])
                else:
                    self.Tw[i, j] = 10

    cpdef double getBiasTrans(self, init, final):
        """
        """
        cdef double delB = final - init
        cdef double prob = getProb(delB, 0, self.sigmaB, 0.00001)

        if prob != 0:
            return log(prob)
        else:
            return -100

    cpdef findSequence(self):
        """
        find most likely sequence of states given observations (obs)
        """

        self.V[0, 0] = math.log(getProb(self.states[0], 0, 0.75, self.resW)
        self.V[0, self.Nw-1] = math.log(getProb(self.states[Nw-1], 0, 0.75, self.resW))

        self.B[0, 0] = 0
        self.B[0, self.Nw-1] = self.Nw-1

        for i in range(1, self.Nw-1):
            self.V[0, i] = math.log(getProb(self.states[i], 0, 0.75, self.resW))
            for j in range(self.Ns):
                self.B[0, i] = i

        cdef double p
        cdef double p_max
        cdef int s_max

        for t in range(1, self.N):  # looping over time
            for i in range(self.Nw):  # looping over possible new states
                p_max = -1e1000
                s_max = -1
                # looping over old states
                for j in range(self.Nw):
                    if self.Tw[j, i] != 10:
                        p = self.V[t-1, j] + self.Tw[j, i]
                        for sens in range(self.Ns):
                            bi = self.obs[sens, t-1] - self.states[j]
                            bf = self.obs[sens, t] - self.states[i]
                            p += self.getBiasTrans(bi, bf)

                        if p > p_max:
                            p_max = p
                            s_max = j

                self.V[t, i] = p_max
                self.B[t, i] = s_max

        cdef int st = np.argmax(self.V[self.N-1, :])

        for t in range(self.N-1, -1, -1):
            self.Bp[t] = self.states[st]
            st = self.B[t, st]

def run(obs):
    cdef Viterbi v
    v = Viterbi.Viterbi(obs)
    v.createMatrices()
    v.findSequence()

    return v.Bp, v.V, v.B, v.states
#    def run(self):
#        """
#        Run the viterbi algorithm, and plot data against GT
#        """
#        self.createMatrices()
#        self.findSequence()

#        plt.plot(self.Bp, 'r')
#        plt.plot(self.omega, 'g')
# if __name__ == "__main__":
# bias, omega, obs = generateData.generate(2, 150000)
