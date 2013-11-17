import numpy as np
import matplotlib.pyplot as plt


def generate(num, samples):
    """
    Generate bias data to be added to actual 'omega'
    """
    sigmaB = [0.0002, 0.0003]
    sigmaW = 0.0085
    #sigmaN = 0.03

    initB = [0.025, -0.040]

    bias = np.empty([num, samples])
    for i in range(num):
        bias[i, 0] = np.random.normal(initB[i], sigmaB[i])

    omega = np.empty(samples)
    omega[0] = 0

    obs = np.empty([num, samples])
    obs[:, 0] = bias[:, 0] + omega[0]

    for i in xrange(1, samples):
        for j in xrange(num):
            bias[j, i] = bias[j, i-1] + np.random.normal(0, sigmaB[j])
        omega[i] = omega[i-1] + np.random.normal(0, sigmaW)

        obs[:, i] = bias[:, i] + omega[i]

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)

if __name__ == "__main__":

    bias, omega, obs = generate(2, 150000)

    plt.figure(1)
    plt.plot(obs[0, :])

    plt.figure(2)
    plt.plot(obs[1, :])

    plt.figure(3)
    plt.plot(np.cumsum(obs[0, :]),'b')
    plt.plot(np.cumsum(obs[1, :]),'r')
    plt.plot(np.cumsum(omega),'g')

    plt.show()
