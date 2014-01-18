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

    def populateInitial(self, omega, obs):

        # omega should be a list/tuple of two values - mean and variance
        self.particles[0] = np.random.normal(omega[0], omega[1], self.Np)

        # bias should be a list of lists/tuples
        for j in xrange(self.Np):
            for i in xrange(self.dim):
                self.particles[i+1, j] =\
                                         np.random.normal(obs[i]-self.particles[0, j],
                                                          self.meas[i][1])

    def draw(self, weights, choices, count):

        N = choices.shape[1]
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

        for i in xrange(self.Np):
            self.particles[0, i] += np.random.normal(self.process[0][0],
                                                     self.process[0][1])

            for j in xrange(self.dim):
                self.particles[j+1, i] = self.process[j+1][0]*self.particles[j+1, i] +\
                                         np.random.normal(0, self.process[j+1][1])

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

    def run(self, obs, omega, b):

        self.setProcessModel((0.0, 0.0085), [(0.8, 0.09)]*self.dim)
        self.setMeasNoise([(0.0, 0.000001)]*self.dim)
        self.populateInitial([omega[0], 0.005], obs[:, 0])

        estimate = np.zeros([obs.shape[1]-1, obs.shape[0]+1])
        Z = obs.T

        N = obs.shape[1]

        for i in xrange(N-1):
            self.predict()
            self.measure(Z[i+1])

            estimate[i] = self.getEstimate()

            self.resample()

        colors = ['r', 'b', 'y', 'c', 'm']

        m = np.copy(obs)
        for i in xrange(m.shape[0]):
            m[i] = m[i] - m[i, 0] + omega[0]

        m = np.mean(m, axis=0)
        MEAN = False

        plt.close()
        plt.figure()

        plt.subplot(211)
        plt.plot(omega[:N], 'g', label="Ground truth")
        plt.plot(estimate[:, 0], colors[0], label="PF Estimate")
        #plt.plot(b.T, label="Bias")
        if MEAN:
            plt.plot(m, 'b', label="Mean")
        plt.title("Particle filter")
        plt.legend()
        plt.xlabel("Time")
        plt.ylabel("Angular velocity")

        plt.subplot(212)
        plt.plot(omega[:N-1]-estimate[:, 0], 'r', label="Estimation Error")
        if MEAN:
            plt.plot(omega[:N-1]-m[:-1], 'b', label="Mean estim. Error")
        plt.title("Particle filter error")
        plt.legend()
        plt.xlabel("Time")
        plt.ylabel("Error in rate estimation")

        plt.show()

    def compare(self, obs, omega):

        plt.close()
        plt.figure(1)

        N = obs.shape[1]

        plt.subplot(211)
        plt.plot(omega[:N], 'g', label="Ground truth")

        for num in xrange(1, obs.shape[0]+1):

            self = ParticleFilter(500, num)

            self.setProcessModel((0.0, 0.0085), [(0.0, 0.00015)]*self.dim)
            self.setMeasNoise([(0.0, 0.005)]*self.dim)
            self.populateInitial([0, 0.1], obs[:num, 0])

            estimate = np.zeros([N-1, self.dim+1])
            Z = obs[:num].T

            for i in xrange(N-1):
                self.predict()
                self.measure(Z[i+1])

                estimate[i] = self.getEstimate()

                self.resample()

            colors = ['r', 'b', 'y', 'c', 'm']

            plt.subplot(211)
            plt.plot(estimate[:, 0], colors[num], label="Estimate %d" % (num+1))

            plt.subplot(212)
            plt.plot(np.abs(omega[:N-1]-estimate[:, 0]), colors[num], label="Error %d" % (num+1))

        plt.subplot(211)
        plt.title("Particle filter")
        plt.legend()
        plt.xlabel("Time")
        plt.ylabel("Angular velocity")

        plt.subplot(212)
        plt.title("Particle filter error")
        plt.legend()
        plt.xlabel("Time")
        plt.ylabel("Error in rate estimation")

        plt.show()
