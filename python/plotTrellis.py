import matplotlib.pyplot as plt
import numpy as np


def find_nearest_idx(states, value):
    idx = (np.abs(states-value)).argmin()
    return idx

def plot(Bp, B, states, omega):

    N = Bp.shape[0]
    minst = np.min(Bp) - 0.005
    maxst = np.max(Bp) + 0.005

    for i in range(1, N):
        ml = find_nearest_idx(states, Bp[i])
        src = B[i, ml]
        source_st = states[src]
        ml_st = Bp[i]

        plt.plot([i-1, i], [source_st, ml_st], 'r')

        minst = min(minst, source_st, ml_st)
        maxst = max(maxst, source_st, ml_st)

    plt.yticks(np.arange(minst, maxst, 0.005))

    plt.show()
