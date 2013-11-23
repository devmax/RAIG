import numpy as np
import matplotlib.pyplot as plt
import math
import itertools


def generate(num, samples):
    """
    Generate bias and omega
    """
    sigmaB = np.random.normal(0.00025, 0.000035, num)
    sigmaW = 0.0085
    #sigmaN = 0.03

    initB = np.random.normal(0.00025, 0.25, num)

    bias = np.zeros([num, samples])
    bias[:, 0] = initB.T

    omega = np.zeros(samples)
    omega[0] = 0

    obs = np.empty([num, samples])

    for i in xrange(num):
        bias[i] = bias[i] + np.cumsum(np.random.normal(0, sigmaB[i], samples))

    obs[:, 0] = bias[:, 0] + omega[0]
    generated = 1

    dw = np.random.normal(0, sigmaW, samples)

    w = np.cumsum(dw)
    idx_t = np.where((-0.01 < w) & (w < 0.01))
    idx = [idx_t[0][10]]
    cur = 2

    while len(idx) < 15:
        if(idx_t[0][cur] - idx[-1] > 1000):
            idx.append(idx_t[0][cur])
        cur += 1
        if cur == idx_t[0].shape[0]:
            break

    cur = 0
    idx = np.concatenate((idx, [samples]))

    while generated < samples:

        if generated < idx[cur]:
            omega[generated] = omega[generated-1] + dw[generated-1]
            obs[:, generated] = bias[:, generated] + omega[generated]
            generated += 1
        else:
            generated = idx[cur+1] 
            obs[:, idx[cur]:generated] = omega[idx[cur]:generated] + bias[:, idx[cur]:generated]
            cur = min(cur + 2, idx.shape[0])

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
