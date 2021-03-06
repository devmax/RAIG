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
    int abs(int n)

cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b

cpdef double logProb(double x, double mu, double sigma):
    return -(log(sigma)+0.5*(((x-mu)/sigma)**2))


@cython.boundscheck(False)
cpdef createMatrices(int Nw, double resW, double sigmaW):

    cdef np.ndarray[np.double_t, ndim = 2] Tw = np.empty([Nw, Nw])

    for j in xrange(Nw):
        for i in xrange(Nw):
            # log-probability that state j came from state i
            Tw[i, j] = logProb(abs(j-i)*resW, 0, sigmaW)

    return Tw

cpdef double getBiasTrans(double init, double final, double sigmaB):

    return logProb(final, 0.8*init, sigmaB)


@cython.boundscheck(False)
cpdef findSequence(double resW, double sigmaW, double sigmaB,
                   np.ndarray[np.double_t, ndim=1] states,
                   np.ndarray[np.double_t, ndim=2] obs):

    cdef int Nw = states.shape[0]
    cdef int N = obs.shape[1]
    cdef int Ns = obs.shape[0]

    cdef np.ndarray[np.double_t, ndim = 2] V = np.ones([2, Nw])
    cdef np.ndarray[np.int_t, ndim = 2] B = np.empty([N, Nw], dtype=np.int)

    cdef np.ndarray[np.double_t, ndim = 1] Bp = np.zeros(N, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim = 2] Tw = createMatrices(Nw,
                                                               resW, sigmaW)

    for i in xrange(Nw):
        V[0, i] = logProb(states[i], 0., 0.5)

    B[0] = np.arange(Nw, dtype=np.int)
    Bp[0] = 0.

    cdef double p
    cdef double p_max
    cdef unsigned int s_max
    cdef double p_ml
    cdef unsigned int ml

    for t in xrange(1, N):  # looping over time
        if t % 1000 == 0:
            print "Observation number: ", t
        p_ml = -1e1000
        for i in xrange(Nw):  # looping over possible new states
            p_max = -1e1000
            # looping over old states
            for j in xrange(Nw):
                p = V[0, < unsigned int > j] + Tw[< unsigned int > j, < unsigned int > i]
                for sens in range(Ns):
                    bi = obs[< unsigned int > sens, < unsigned int
                             > (t-1)] - states[< unsigned int >
                                               j]
                    bf = obs[< unsigned int > sens, < unsigned int > t] - states[< unsigned int > i]
                    p += getBiasTrans(bi, bf, sigmaB)

                if p > p_max:
                    p_max = p
                    s_max = j

            V[1, < unsigned int > i] = p_max
            B[< unsigned int > t, < unsigned int > i] = s_max

            if p_max > p_ml:
                p_ml = p_max
                ml = i
                Bp[< unsigned int > t] = states[< unsigned int > ml]

        V[0] = np.copy(V[1])

    cdef np.ndarray[np.double_t, ndim = 1] Bpv = np.zeros(N, dtype=np.double)

    for t in xrange(N-1, -1, -1):
        Bpv[< unsigned int > t] = states[< unsigned int > ml]
        ml = B[< unsigned int > t, < unsigned int > ml]

    return Bp, Bpv, V, B


def estimate(obs, omega, resW=0.005, sigmaB=0.09):

    N = obs.shape[1]
    omega = omega[:N]

    sigmaW = 0.0085
    print "bias standard deviation = ", sigmaB
    print "Resolution = ", resW

    Ns = obs.shape[0]

    lim = double_max((np.min(obs)*-1), np.max(obs))
    states = np.concatenate((np.arange(0, -lim, -resW)[:0:-1],
                             np.arange(0, lim, resW)), 1)

    print states.shape[0], "states, [", states[0], ",", states[-1], "]"

    Bp, Bpv, V, B = findSequence(resW, sigmaW, sigmaB, states, obs)

    plt.figure(1)
    plt.subplot(211)

    plt.plot(Bpv, color='r', label="Estimated rate(backtrack)")
    plt.plot(omega, color='g', label="True angular rate")
    plt.legend()
    plt.title("Viterbi estimate")
    plt.xlabel('Time')
    plt.ylabel('Angular Rate')

    plt.subplot(212)

    plt.plot(Bp, color='r', label="Estimated rate(forward)")
    plt.plot(omega, color='g', label="True angular rate")
    plt.title("Instantaneous estimates")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Angular Rate")

    plt.figure(2)
    plt.subplot(211)

    err = Bp - omega

    plt.plot(err, label="Error in angular rate estimation")
    #plt.plot(np.cumsum(err), label="Cumulative error in rate")
    plt.legend()
    plt.title("Instantaneous estimate error")
    plt.xlabel('Time')
    plt.ylabel('Error')

    plt.subplot(212)
    err = Bpv - omega

    plt.plot(err, label="Error in angular rate estimation")
    #plt.plot(np.cumsum(err), label="Cumulative error in rate")
    plt.legend()
    plt.title("Viterbi estimate error")
    plt.xlabel('Time')
    plt.ylabel('Error')

    plt.show()

    return Bp, Bpv, V, B, states, err


def run(obs, omega):

    print "Using true sigma values"
    N = obs.shape[1]
    omega = omega[:N]

    count = 0
    res = []
    numStates = []
    error = []
    estimates = []
    resW = 1.

    res = [1.0, 0.5, 0.25, 0.125, 0.06, 0.04, 0.02, 0.01, 0.008, 0.005]#,
           #0.0025, 0.0012, 0.00075, 0.0005, 0.0002, 0.0001, 0.00005]

    sigmaW = 0.0085
    sigmaB = 0.00015

    Ns = obs.shape[0]

    for resW in res:
        print "\nWorking on resolution of ", resW
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

        Bp, Bpf, V, B = findSequence(resW, sigmaW, sigmaB, states, obs)
        err = Bpf - omega

        numStates.append(states.shape[0])
        error.append(err)
        # estimates.append(Bpf)
        count += 1
        print "\n", states.shape[0], " states, error= ",\
            np.sqrt(np.sum(np.power(err, 2))/N)
        print "*"*25

    return numStates, res, error, estimates
