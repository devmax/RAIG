import numpy as np
import matplotlib.pyplot as plt


def generateVariant(num, samples, lim):
    """
    Generate bias and omega, but with sharper penalties for omega at the limits
    """
    sigmaB = np.random.normal(0.00015, 0.000015, num)
    sigmaW = 0.0085

    initB = np.random.normal(0.00025, 0.025, num)

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

        #sig = (1 - (omega[generated-1]/lim))*sigmaW
        sig = sigmaW
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
            bias[:, generated] = bias[:, generated-1] + dbias[:,
                                                              generated-1]
            obs[:, generated] = bias[:, generated] + omega[generated]
            generated += 1

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)


def sanityCheck(num, samples):
    """
    Generate bias and omega, with dw being quantized by to factors of
    0.005

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

    dw = np.random.normal(0, sigmaW, samples)

    while generated < samples:
        delta = round(dw[generated-1]*1000.)
        r = delta % 5
        delta += 5-r if r > 2 else -r

        omega[generated] = omega[generated-1] + delta/1000.
        bias[:, generated] = bias[:, generated-1] + dbias[:, generated-1]
        obs[:, generated] = bias[:, generated] + omega[generated]
        generated += 1

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)

if __name__ == "__main__":

    num = 4

    bias, omega, obs = generateVariant(num, 500000, 1.4)
    #bias, omega, obs = sanityCheck(num, 150000)
    #bias, omega, obs = generate(num, 150000)

    plt.figure()
    plt.plot(obs.T)
    plt.title('Observations')
    plt.xlabel('Time')
    plt.ylabel('Observations')

    plt.figure()
    plt.plot(np.cumsum(omega), 'g')
    plt.plot(np.cumsum(obs, axis=1).T)
    plt.ylabel('Angle')
    plt.title('Integrated angular rate')

    plt.figure()
    plt.plot(omega, 'g')
    plt.ylabel('True omega')
    plt.title('Ground truth')

    plt.show()
