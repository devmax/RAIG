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
    double round(double x)

cdef inline double double_max(double a, double b): return a if a>=b else b
cdef inline double double_min(double a, double b): return a if a<=b else b

cpdef double getPDF(double sqrt_tpi, double x, double mu, double sigma):
    return (1/(sqrt_tpi*sigma))*exp((-0.5)*pow(((x-mu)/sigma), 2))

cpdef double getProb(double x, double mu, double sigma, double res):
    return (0.5*(1+erf((x-mu+res)/(sqrt(2)*sigma))) -
            0.5*(1+erf((x-mu-res)/(sqrt(2)*sigma))))


@cython.boundscheck(False)
cpdef createMatrices(int Nw, double resW, double sigmaW):

    cdef double sqrt_tpi = sqrt(2*math.pi)

    cdef np.ndarray[np.double_t, ndim=2] Tw = np.empty([Nw, Nw])

    for i in range(Nw):
        for j in range(Nw):
            Tw[<unsigned int>i, <unsigned int>j] = getPDF(sqrt_tpi,
                                                          (j-i)*resW,
                                                          0, sigmaW)/100.0

        for j in range(Nw):
            if Tw[<unsigned int>i, <unsigned int>j] != 0:
                Tw[<unsigned int>i, <unsigned int>j] =
                log(Tw[<unsigned int>i, <unsigned int>j])
            else:
                Tw[<unsigned int>i, <unsigned int>j] = 1.0

    return Tw

cpdef double getBiasTrans(double init, double final, double sigmaB):

    cdef double delB = final - init
    cdef double prob = getProb(delB, 0, sigmaB, 0.00001)

    if prob != 0:
        return log(prob)
    else:
        return -9999


@cython.boundscheck(False)
cpdef findSequence(double resW, double sigmaW, double sigmaB,
                   np.ndarray[np.double_t, ndim=1] states,
                   np.ndarray[np.double_t, ndim=2] obs):

    cdef int Nw = states.shape[0]
    cdef int N = obs.shape[1]
    cdef int Ns = obs.shape[0]

    cdef np.ndarray[np.double_t, ndim = 2] V = np.ones([N, Nw])*(10.0)
    cdef np.ndarray[np.int_t, ndim = 2] B = np.empty([N, Nw], dtype=np.int)

    cdef np.ndarray[np.double_t, ndim = 1] Bp = np.zeros(N, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim = 2] Tw = createMatrices(Nw, resW, sigmaW)

    cdef sqrt_tpi = sqrt(2*math.pi)

    for i in xrange(Nw):
        V[0, i] = getProb(states[i], 0, 0.5, resW)

    B[0] = np.arange(Nw, dtype=np.int)
    Bp[0] = 0.

    cdef double p
    cdef double p_max
    cdef unsigned int s_max
    cdef double p_ml
    cdef unsigned int ml

    cdef double oi
    cdef double of
    cdef float r
    cdef double div = 1000.0

    for t in xrange(1, N):  # looping over time
        p_ml = -1e1000
        for i in xrange(Nw):  # looping over possible new states
            p_max = -1e1000
            # looping over old states
            for j in xrange(Nw):
                if Tw[<unsigned int> j, <unsigned int> i] <= 0.:
                    p = V[<unsigned int> (t-1), <unsigned int> j] + Tw[<unsigned int> j, <unsigned int> i]
                    for sens in range(Ns):
                        oi = round(obs[< unsigned int > sens, < unsigned
                                       int > (t-1)] * div)
                        of = round(obs[<unsigned int> sens, <unsigned
                                       int> t] * div)
                        r = oi % 5
                        oi += 5-r if r > 2 else -r

                        r = of % 5
                        of += 5-r if r > 2 else -r

                        bi = (oi/div) - states[<unsigned int>j]
                        bf = (of/div) - states[<unsigned int>i]
                        p += getBiasTrans(bi, bf, sigmaB)

                    if p > p_max:
                        p_max = p
                        s_max = j

            V[<unsigned int>t, <unsigned int>i] = p_max
            B[<unsigned int>t, <unsigned int>i] = s_max

            if p_max > p_ml:
                p_ml = p_max
                ml = i
        Bp[<unsigned int>t] = states[<unsigned int>ml]

    #cdef double vmax = -1e100000000
    #cdef unsigned int st

    #for i in xrange(Nw):
    #    if(V[<unsigned int>(N-1), <unsigned int>i] < 0
    #       and V[<unsigned int>(N-1), <unsigned int>i] > vmax):
    #        vmax = V[<unsigned int>(N-1), <unsigned int>i]
    #        st = i

    #for t in xrange(N-1, -1, -1):
    #    Bp[<unsigned int>t] = states[<unsigned int>st]
    #    st = B[<unsigned int>t, <unsigned int>st]

    return Bp, V, B


def run(obs):

    print "Just checking"
    resW = 0.005
    sigmaW = 0.0085
    sigmaB = 0.0025

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

    states = np.concatenate((np.arange(0, -lim, -resW)[:0:-1],
                             np.arange(0, lim, resW)), 1)

    Bp, V, B = findSequence(resW, sigmaW, sigmaB, states, obs)

    return Bp, V, B, states
