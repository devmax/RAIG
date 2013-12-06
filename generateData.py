import numpy as np
import matplotlib.pyplot as plt


def generateVariant(num, samples, lim):
    """
    Generate bias and omega, but with sharper penalties for omega at the limits
    """
    sigmaB = np.random.normal(0.00015, 0.000015, num)
    sigmaW = 0.0085

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

    while generated < samples:

        sig = (1 - (omega[generated-1]/lim))*sigmaW
        # mu = -omega[generated-1]
        dw = np.random.normal(0, sig, 1)
        omega[generated] = omega[generated-1] + dw
        if abs(omega[generated]) > lim:
            omega[generated] -= 2*dw
        bias[:, generated] = bias[:, generated-1] + dbias[:, generated-1]
        obs[:, generated] = bias[:, generated] + omega[generated]
        generated += 1

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)


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

    bias, omega, obs = generateVariant(num, 150000, 1.8)

    for i in xrange(num):
        plt.figure(i+1)
        plt.plot(obs[i, :])

        plt.figure(num+1)
        plt.plot(np.cumsum(obs[i, :]))

    plt.figure(num+1)
    plt.plot(np.cumsum(omega), 'g')

    plt.figure(10)
    plt.plot(omega, 'r')

    plt.show()
