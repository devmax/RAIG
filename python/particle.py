import numpy as np
import math
import matplotlib.pyplot as plt


class ParticleFilter:

    def __init__(self, Np, dim):

        self.Np = Np
        self.dim = dim

        self.particles = np.empty([self.dim+1, self.Np], dtype=float)
        self.w = np.ones(self.Np)

        self.process = []

    def setProcessModel(self, omega, bias):

        self.process = [omega]  # omega is a tuple
        self.process.extend(bias)  # bias is a 'dim'-d list of tuples

    def setMeasNoise(self, meas):

        self.meas = meas  # meas should be a 'dim'-d list of tuples

    def populateInitial(self, omega, bias):

        # omega should be a list/tuple of two values - mean and variance
        self.particles[0] = np.random.normal(omega[0], omega[1], self.Np)

        # bias should be a list of lists/tuples
        for i in xrange(self.dim):
            self.particles[(i+1)] = np.random.normal(bias[i][0],
                                                     bias[i][1], self.Np)

    def draw(self, weights, choices, count):

        N = choices.shape[0]
        w_max = max(weights)
        idx = np.random.randint(0, N-1)
        beta = 0

        samples = np.empty_like(choices)

        for i in xrange(count):
            beta += np.random.uniform(0, 2*w_max)
            while weights[idx] < beta:
                beta -= weights[idx]
                idx = (idx + 1) % N

            samples[:, i] = choices[:, idx]

        return samples

    def predict(self):

        self.particles[0] += np.random.normal(self.process[0][0],
                                              self.process[0][1], self.Np)

        for i in xrange(self.dim):
            self.particles[i+1] += np.random.normal(self.process[i+1][0],
                                                    self.process[i+1][1],
                                                    self.Np)

    def Gaussian(self, x, mu, sigma):
        return (1/(math.sqrt(2)*sigma))*math.exp((-0.5)*pow(((x-mu)/sigma), 2))

    def measure(self, Z):

        assert Z.shape[0] == self.dim

        for j in xrange(self.Np):
            for i in xrange(self.dim):
                x = Z[i]-self.particles[0, j]-self.particles[i+1, j]
                p = self.Gaussian(x, self.meas[i][0], self.meas[i][1])
                self.w[j] = self.w[j]*p

    def getEstimate(self):

        return self.particles[:, np.argmax(self.w)]

    def resample(self):

        self.particles = np.copy(self.draw(self.w, self.particles, self.Np))
        self.w = np.ones(self.Np)*1./self.Np*1.


def runFilter(obs, omega):

    pf = ParticleFilter(50, obs.shape[0])

    pf.populateInitial([omega[0], 0.05], [[obs[i, 0]-omega[0], 0.01]
                                          for i in xrange(pf.dim)])
    pf.setProcessModel((0.0, 0.0085), [(0.0, 0.00015)]*pf.dim)
    pf.setMeasNoise([(0.0, 0.001)]*pf.dim)

    estimate = np.empty([obs.shape[1], obs.shape[0]+1])
    Z = obs.T

    N = obs.shape[1]

    for i in xrange(N):
        pf.predict()
        pf.measure(Z[i])

        estimate[i] = pf.getEstimate()

        pf.resample()

    plt.plot(omega[:N], 'g')
    plt.plot(estimate[:, 0], 'r')

    plt.show()
