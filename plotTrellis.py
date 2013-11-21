import matplotlib.pyplot as plt
import numpy as np


def plot(V, B, states, omega):
    N = V.shape[0]
    Nw = V.shape[1]

    Bp = np.empty(N)

    for i in xrange(N):
        max_val = -1e100000000
        max_idx = None
        for j in xrange(Nw):
            if (V[i, j] < 0 and V[i, j] > max_val):
                max_val = V[i, j]
                max_idx = j

        Bp[i] = states[max_idx]

    plt.plot(Bp, 'r')
    plt.plot(omega[:N], 'g')

    plt.show()
