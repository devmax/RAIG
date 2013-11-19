import numpy as np
import scipy.stats as stats
import math
cimport numpy as np
cimport cython

np.import_array()

import matplotlib.pyplot as plt

cdef extern from "math.h":
    double exp(double x)
    double log(double x)
    double erf(double x)
    double sqrt(double x)
    double pow(double base, double exponent)

cdef inline double double_max(double a, double b): return a if a>=b else b
cdef inline double double_min(double a, double b): return a if a<=b else b

cpdef double getPDF(double sqrt_tpi, double x, double mu, double sigma):
    return (1/(sqrt_tpi*sigma))*exp((-0.5)*pow(((x-mu)/sigma), 2))

cpdef double getProb(double x, double mu, double sigma, double res):
    return (0.5*(1+erf((x+res-mu)/(sqrt(2)*sigma))) -
            0.5*(1+erf((x-res-mu)/(sqrt(2)*sigma))))

cpdef createMatrices(int Nw, double resW, double sigmaW):

    cdef double sqrt_tpi = sqrt(2*math.pi)

    cdef np.ndarray[np.double_t, ndim=2] Tw = np.empty([Nw, Nw])

    for i in range(Nw):
        for j in range(Nw):
            Tw[i, j] = getPDF(sqrt_tpi, (j-i)*resW, 0, sigmaW)
        Tw[i, :] /= np.sum(Tw[i, :])
        for j in range(Nw):
            if Tw[i, j] != 0:
                Tw[i, j] = log(Tw[i, j])
            else:
                Tw[i, j] = 10.0

    return Tw

cpdef double getBiasTrans(double init,double final, double sigmaB):

    cdef double delB = final - init
    cdef double prob = getProb(delB, 0, sigmaB, 0.00001)

    if prob != 0:
        return log(prob)
    else:
        return -100

cpdef findSequence(double resW, double sigmaW, double sigmaB,
                   np.ndarray[np.double_t, ndim=1] states, np.ndarray[np.double_t, ndim=2] obs):

    cdef int Nw = states.shape[0]
    cdef int N = obs.shape[1]
    cdef int Ns = obs.shape[0]

    cdef np.ndarray[np.double_t, ndim=2] V = np.ones([N, Nw])*-1500
    cdef np.ndarray[np.int_t, ndim=2] B = np.empty([N, Nw], dtype=np.int)

    cdef np.ndarray[np.int_t, ndim=1] Bp = np.empty(N, dtype=np.int)

    cdef np.ndarray[np.double_t, ndim=2] Tw = createMatrices(Nw, resW, sigmaW)

    cdef double prob = getProb(states[0], 0, 0.75, resW)
    if prob != 0:
        V[0, 0] = math.log(prob)
    else:
        V[0, 0] = -150

    prob = getProb(states[Nw-1], 0, 0.75, resW)
    if prob != 0:
        V[0, Nw-1] = math.log(prob)
    else:
        V[0, Nw-1] = -150

    B[0, 0] = 0
    B[0, Nw-1] = Nw-1

    for i in range(1, Nw-1):
        prob = getProb(states[i], 0, 0.75, resW)
        if prob != 0:
            V[0, i] = math.log(prob)
        else:
            V[0, i] = -150

        for j in range(Ns):
            B[0, i] = i

    cdef double p
    cdef double p_max
    cdef int s_max

    for t in range(1, N):  # looping over time
        for i in range(Nw):  # looping over possible new states
            p_max = -1e1000
            s_max = -1
            # looping over old states
            for j in range(Nw):
                if Tw[j, i] != 10:
                    p = V[t-1, j] + Tw[j, i]
                    for sens in range(Ns):
                        bi = obs[sens, t-1] - states[j]
                        bf = obs[sens, t] - states[i]
                        p += getBiasTrans(bi, bf, sigmaB)

                    if p > p_max:
                        p_max = p
                        s_max = j

            V[t, i] = p_max
            B[t, i] = s_max

    cdef int st = np.argmax(V[N-1, :])

    for t in range(N-1, -1, -1):
        Bp[t] = states[st]
        st = B[t, st]

    return Bp, V, B

def run(obs):

    resW = 0.005
    sigmaW = 0.0085
    sigmaB = 0.00025

    Ns = obs.shape[0]

    minw = 1e10
    maxw = -1e10

    for i in xrange(Ns):
        val = np.amin(obs[i, :])
        if val < minw:
            minw = val

        val = np.amax(obs[i, :])
        if val > maxw:
            maxw = val

    lim = double_max((minw*-1), maxw)

    states = np.concatenate((np.arange(-lim, 0, resW),
                                      np.arange(0, lim, resW)), 1)

    Bp, V, B = findSequence(resW, sigmaW, sigmaB, states, obs)

    return Bp, V, B, states
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
