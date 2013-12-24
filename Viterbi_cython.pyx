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

cpdef double getPDF(double sqrt_tpi, double x, double mu, double sigma):
    return (1/(sqrt_tpi*sigma))*exp((-0.5)*pow(((x-mu)/sigma), 2))

cpdef double getProb(double x, double mu, double sigma, double res):
    return (0.5*(1+erf((x-mu+res)/(sqrt(2)*sigma))) -
            0.5*(1+erf((x-mu-res)/(sqrt(2)*sigma))))


@cython.boundscheck(False)
cpdef createMatrices(int Nw, double resW, double sigmaW):

    cdef double sqrt_tpi = sqrt(2*math.pi)

    cdef np.ndarray[np.double_t, ndim = 2] Tw = np.empty([Nw, Nw])

    for j in xrange(Nw):
        for i in xrange(Nw):
            # probability that state j came from state i
            Tw[<unsigned int>i, <unsigned int>j] = getPDF(sqrt_tpi,
                                                          abs(j-i)*resW,
                                                          0, sigmaW)
        Tw[:, j] /= 100. # sum(Tw[:, j])
        for i in range(Nw):
            p = Tw[<unsigned int>i, <unsigned int>j]
            if p != 0:
                Tw[<unsigned int>i, <unsigned int>j] = log(p)
            else:
                Tw[<unsigned int>i, <unsigned int>j] = 10.0
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

    cdef np.ndarray[np.double_t, ndim = 2] V = np.ones([2, Nw])*(10.0)
    cdef np.ndarray[np.int_t, ndim = 2] B = np.empty([N, Nw], dtype=np.int)

    cdef np.ndarray[np.double_t, ndim = 1] Bp = np.zeros(N, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim = 2] Tw = createMatrices(Nw,
                                                               resW, sigmaW)

    cdef double sqrt_tpi = sqrt(2*math.pi)
    cdef double pr

    for i in xrange(Nw):
        pr = getProb(states[i], 0, 0.5, resW)
        if pr != 0:
            V[0, i] = log(pr)
        else:
            V[0, i] = 1e2

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
                if Tw[< unsigned int > j, < unsigned int > i] <=\
                0. and V[0, < unsigned int > j] <= 0:
                    p = V[0, < unsigned int > j]\
                    + Tw[< unsigned int > j, < unsigned int > i]
                    for sens in range(Ns):
                        bi = obs[< unsigned int > sens, < unsigned int
                                 > (t-1)] - (states[< unsigned int >
                                                    j])
                        bf = obs[< unsigned int > sens, < unsigned int
                                 > t] - (states[< unsigned int > i])
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


def estimate(obs, omega, resW, sigmaB):

    print "Randomly offsetting states, foo!"

    N = obs.shape[1]
    omega = omega[:N]

    sigmaW = 0.0085
    print "bias standard deviation = ", sigmaB
    print "Resolution = ", resW

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

    print states.shape[0], "states, [", states[0], ",", states[-1], "]"

    Bp, Bpv, V, B = findSequence(resW, sigmaW, sigmaB, states, obs)

    plt.figure(1)
    plt.subplot(211)

    plt.plot(Bpv, color='r', label="Estimated rate(backtrack)")
    plt.plot(omega, color='g', label="True angular rate")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2,
               mode="expand", borderaxespad=0.)
    plt.title("Viterbi estimate")
    plt.xlabel('Time')
    plt.ylabel('Angular Rate')

    plt.subplot(212)

    plt.plot(Bp, color='r', label="Estimated rate(forward)")
    plt.plot(omega, color='g', label="True angular rate")
    plt.title("Instantaneous estimates")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2,
               mode="expand", borderaxespad=0.)
    plt.xlabel("Time")
    plt.ylabel("Angular Rate")

    plt.figure(2)

    plt.subplot(211)
    err = Bp - omega

    plt.plot(err, label="Error in angular rate estimation")
    #plt.plot(np.cumsum(err), label="Cumulative error in rate")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2,
               mode="expand", borderaxespad=0.)
    plt.title("Instantaneous estimate error")
    plt.xlabel('Time')
    plt.ylabel('Error')

    plt.subplot(212)
    err = Bpv - omega

    plt.plot(err, label="Error in angular rate estimation")
    #plt.plot(np.cumsum(err), label="Cumulative error in rate")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2,
               mode="expand", borderaxespad=0.)
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
