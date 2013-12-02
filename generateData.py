import numpy as np
import matplotlib.pyplot as plt


def generate(num, samples):
    """
    Generate bias and omega
    """
    sigmaB = np.random.normal(0.00015, 0.000015, num)
    sigmaW = 0.0085
    #sigmaN = 0.03

    initB = np.random.normal(0.00025, 0.25, num)

    bias = np.zeros([num, samples])
    bias[:, 0] = initB.T
    dbias = np.zeros([num, samples])

    omega = np.zeros(samples)
    omega[0] = 0

    obs = np.empty([num, samples])

    for i in xrange(num):
        dbias[i] = np.random.normal(0, sigmaB[i], samples)

    obs[:, 0] = bias[:, 0] + omega[0]
    generated = 1

    dw = np.random.normal(0, sigmaW, samples)

    while generated < samples:

            omega[generated] = omega[generated-1] + dw[generated-1]
            bias[:, generated] = bias[:, generated-1] + dbias[:, generated-1]
            obs[:, generated] = bias[:, generated] + omega[generated]
            generated += 1

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)

if __name__ == "__main__":

    num = 4

    bias, omega, obs = generate(num, 150000)

    for i in xrange(num):
        plt.figure(i+1)
        plt.plot(obs[i, :])

        plt.figure(num+1)
        plt.plot(np.cumsum(obs[i, :]))

    plt.figure(num+1)
    plt.plot(np.cumsum(omega), 'g')

    plt.show()
