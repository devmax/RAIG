import numpy as np
import math


class ParticleFilter:

    particles = None
    w = None

    def __init__(self, N, dim):

        self.N = N
        self.dim = dim

        particles = np.empty([self.dim, self.N], dtype=float)
        w = np.empty(self.N, dtype=float)

    def setProcessModel(self, omega, bias):

        self.process = [omega]  # omega is a tuple
        self.process.extend(bias)  # bias is a 'dim'-d list of tuples

    def setMeasNoise(self, meas):

        self.meas = meas  # meas should be a 'dim'-d list of tuples

    def populateInitial(self, omega, bias):

        self.particles[0] = np.random.normal(omega[0], omega[1], self.N)

        for i in xrange(self.dim-1):
            self.particles[(i+1)] = np.random.normal(bias[i][0],
                                                     bias[i][1], self.N)

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

            samples[i] = choices[idx]

        return samples

    def predict(self):

        self.particles[0] += np.random.normal(self.process[0][0],
                                              self.process[0][1], self.N)

        for i in xrange(self.dim-1):
            self.particles[i+1] += np.random.normal(self.bias[i+1][0],
                                                    self.bias[i+1][1], self.N)

    def Gaussian(self, x, mu, sigma):
        return math.exp((-0.5)*pow(((x-mu)/sigma),
                                   2))/(math.sqrt(2)*sigma*np.pi)

    def measure(self, Z):

        assert Z.dim == self.dim-1

        for i in xrange(self.N):
            self.w[i] = self.Gaussian(Z[0]-self.particles[0] -
                                      self.particles[i+1], self.meas[i][0],
                                      self.meas[i][1])

    def getEstimate(self):

        return self.particles[np.argmax(self.w)]

    def resample(self):

        self.particles = np.copy(self.draw(self.w, self.particles, self.N))


def runFilter(obs, omega):

    pf = ParticleFilter(obs.shape[1], obs.shape[0])

    pf.setProcessModel((0.0, 0.0085), [(0.0, 0.00015)]*pf.dim)
    pf.setMeasNoise((0.0, 0.00001)*pf.dim)

    
