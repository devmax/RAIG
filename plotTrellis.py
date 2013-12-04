import matplotlib.pyplot as plt
import numpy as np


def plot(V, B, states, omega):

    N = V.shape[0]

    for i in range(1, N):
        ml = np.argmax(V[i])
        src = B[i, ml]

        plt.plot([i-1, states[src]], [i, states[ml]], 'r')

    plt.set_yticks(states)
