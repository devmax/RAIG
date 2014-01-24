import numpy as np
import math
import matplotlib.pyplot as plt


class Particle:

    def __init__(self, w, b):

        self.dim = b.shape[0]

        self.w = w
        self.b = b

        self.prev = [[0.0, 0.0, 0.0, 0.0, 0.0]]*self.dim

    def predict(self, dw, db):

        self.w += dw
        self.b += db

        for i in xrange(self.dim):
            self.prev[i].pop(0)
            self.prev[i].append(db[i])

    def getMean(self):

        return np.array([np.mean(self.prev[i]) for i in xrange(self.dim)])

    def getHyp(self):

        return self.w + self.b


class ParticleFilter:

    def __init__(self, Np, dim):

        self.Np = Np
        self.dim = dim

        self.particles = []

        self.w = np.ones(self.Np)

        self.process = []

    def setProcessModel(self, omega, bias):

        self.process = [omega]  # omega is a tuple
        self.process.append(bias)  # bias is a 'dim'-d list of tuples

    def setMeasNoise(self, meas):

        self.meas = meas  # meas should be a float, which is the
                          # variance of the noise Gaussian

    def populateInitial(self, omega, obs):

        self.particles = []

        for i in xrange(self.Np):
            w = np.random.normal(omega[0], omega[1])
            b = np.random.normal(obs-w, self.meas)

            self.particles.append(Particle(w, b))

    def draw(self, weights, choices, count):

        N = len(choices)
        w_max = max(weights)
        idx = np.random.randint(0, N-1)
        beta = 0

        samples = []

        for i in xrange(count):
            beta += np.random.uniform(0, 2*w_max)
            while weights[idx] < beta:
                beta -= weights[idx]
                idx = (idx + 1) % N

            samples.append(choices[idx])

        return samples

    def predict(self):

        for i in xrange(self.Np):

            dw = np.random.normal(self.process[0][0], self.process[0][1])
            db = self.process[1][0]*self.particles[i].b
            db += self.process[1][1]*self.particles[i].getMean()

            self.particles[i].predict(dw, db)

    def Gaussian(self, x, mu, sigma):
        try:
            return (1/(math.sqrt(2)*sigma))*math.exp((-0.5)*pow(((x-mu)/sigma),
                                                                2))
        except OverflowError:
            print "(", x, ";", mu, ",", sigma, ")"

    def measure(self, Z):

        assert Z.shape[0] == self.dim

        g = np.vectorize(self.Gaussian)

        for j in xrange(self.Np):
            p = g(Z, self.particles[j].getHyp(), self.meas)

            self.w[j] *= np.prod(p)

    def getEstimate(self):

        idx = np.argmax(self.w)

        return self.particles[idx].w

    def resample(self):

        self.particles = self.draw(self.w, self.particles, self.Np)
        self.w = np.ones(self.Np)

    def run(self, obs, omega):

        self.setProcessModel((0.0, 0.0085), (-0.003, 0.0, 0.0067))
        self.setMeasNoise(0.095)
        self.populateInitial([omega[0], 0.005], obs[:, 0])

        estimate = np.zeros([obs.shape[1]-1])
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
        plt.plot(estimate, colors[0], label="PF Estimate")
        #plt.plot(b.T, label="Bias")
        if MEAN:
            plt.plot(m, 'b', label="Mean")
        plt.title("Particle filter")
        plt.legend()
        plt.xlabel("Time")
        plt.ylabel("Angular velocity")

        plt.subplot(212)
        plt.plot(omega[:N-1]-estimate, 'r', label="Estimation Error")
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

            self.setProcessModel((0.0, 0.0085), [(-0.003, 0.0085)]*self.dim)
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
