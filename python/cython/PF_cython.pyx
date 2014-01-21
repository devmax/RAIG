import numpy as np
cimport numpy as np
cimport cython

np.import_array()

import matplotlib.pyplot as plt

cdef extern from "math.h":
    double exp(double x)
    double log(double x)
    double sqrt(double x)
    double pow(double base, double exp)

cpdef double Gaussian(double x, double mu, double sigma):
    return (1/(sqrt(2.0)*sigma))*exp((-0.5)*pow(((x-mu)/sigma), 2))

@cython.boundscheck(False)
cpdef populateInitial(int Np, int dim, double mu0, double sig0, np.ndarray[np.double_t, ndim=1] Z):

    cdef np.ndarray[np.double_t, ndim = 2] particles = np.zeros([dim+1, Np], dtype=np.double)

    particles[0] = np.random.normal(mu0, sig0, Np)

    for j in xrange(Np):
        for i in xrange(dim):
            particles[i+1, j] = np.random.normal(Z[i]-particles[0, j], sig0, 1)

    return particles

@cython.boundscheck(False)
cpdef draw(np.ndarray[np.double_t, ndim=1] w, np.ndarray[np.double_t, ndim=2] choices):

    cdef double w_max = max(w)
    cdef int N = choices.shape[1]
    idx = np.random.randint(0, N-1)
    cdef double beta = 0

    cdef np.ndarray[np.double_t, ndim=2] particles = np.empty_like(choices)

    for i in xrange(N):
        beta += np.random.uniform(0, 2*w_max)
        while w[idx] < beta:
            beta -= w[idx]
            idx = (idx+1) % N

        particles[:, i] = choices[:, idx]

    return particles

@cython.boundscheck(False)
cpdef run(int Np, np.ndarray[np.double_t, ndim=2] obs,
          np.ndarray[np.double_t, ndim=1] omega, double sig0):

    cdef double rsig = 0.0085
    cdef double biasmu = 0.8
    cdef double biassig = 0.09
    cdef double meassig = 0.009

    cdef int dim = obs.shape[0]
    cdef int N = obs.shape[1]

    cdef double mu0 = omega[0]

    cdef np.ndarray[np.double_t, ndim = 2] particles = populateInitial(Np, dim, mu0, sig0, obs[:, 0])
    cdef np.ndarray[np.double_t, ndim = 1] w = np.ones(Np, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim = 1] estimate = np.zeros(N, dtype=np.double)

    for i in xrange(N-1):

        for j in xrange(Np):
            # predict
            particles[0, j] += np.random.normal(0, rsig)
            for k in xrange(dim):
                particles[k+1, j] *= biasmu
                particles[k+1, j] += np.random.normal(0, biassig)

                #measure
                w[j] *= Gaussian(obs[k, i+1], particles[0, j]+particles[k+1, j], meassig)

        estimate[i+1] = particles[0, np.argmax(w)]

        particles = draw(w, particles)


    plt.close()
    plt.figure()

    plt.subplot(211)
    plt.plot(omega[:N], 'g', label="Ground truth")
    plt.plot(estimate, 'r', label="PF Estimate")
    plt.title("Particle filter")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Angular velocity")

    plt.subplot(212)
    plt.plot(omega[:N]-estimate, 'r', label="Estimation Error")
    plt.title("Particle filter error")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Error in rate estimation")

    plt.show()
