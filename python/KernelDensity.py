import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt


def generatePDF(data):
    dens = sm.nonparametric.KDEMultivariate(data=data.T, var_type='c',
                                            bw='normal_reference')

    points = np.arange(-1.8, 1.8, 0.005)
    pdf = dens.pdf(points)

    plt.plot(points, pdf)
    plt.xlabel('Angular rate')
    plt.ylabel('P(w)')

    plt.show()
