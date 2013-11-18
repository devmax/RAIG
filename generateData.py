import numpy as np
import matplotlib.pyplot as plt


def generate(num, samples):
    """
    Generate bias data to be added to actual 'omega'
    """
    clamp = lambda n, low, high: max(min(high, n), low)

    sigmaB = [0.0002, 0.0003]
    sigmaW = 0.0085
    #sigmaN = 0.03

    max_dw = 0.035
    min_dw = 0.005

    max_db = 0.001
    min_db = 0.00001

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
            db = np.random.normal(0, sigmaB[j])
            db = clamp(db, -max_db, max_db)
            if abs(db) < min_db:
                db = 0.
            bias[j, i] = bias[j, i-1] + db

        dw = np.random.normal(0, sigmaW)
        dw = clamp(dw, -max_dw, max_dw)
        if abs(dw) < min_dw:
            dw = 0.
        omega[i] = omega[i-1] + dw

        obs[:, i] = bias[:, i] + omega[i]

    return bias, omega, obs  # , np.random.normal(0, sigmaN, samples)

if __name__ == "__main__":

    bias, omega, obs = generate(2, 150000)

    plt.figure(1)
    plt.plot(obs[0, :])

    plt.figure(2)
    plt.plot(obs[1, :])

    plt.figure(3)
    plt.plot(np.cumsum(obs[0, :]), 'b')
    plt.plot(np.cumsum(obs[1, :]), 'r')
    plt.plot(np.cumsum(omega), 'g')

    plt.show()
