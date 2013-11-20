import matplotlib.pyplot as plt


def plot(B):

    plt.figure(1)

    for i in xrange(1, B.shape[0]):
        for j in xrange(B.shape[1]):
            plt.plot([i-1, i], [B[i, j], j], 'r')

    plt.show()
