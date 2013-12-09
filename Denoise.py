from statsmodels.robust import stand_mad
import statsmodels.api as sm
import pywt
import numpy as np
import matplotlib.pyplot as plt
import math


def coef_pyramid_plot(coefs, first=0, scale='uniform', ax=None):
    """
    Parameters
    ----------
    coefs : array-like
        Wavelet Coefficients. Expects an iterable in order Cdn, Cdn-1, ...,
        Cd1, Cd0.
    first : int, optional
        The first level to plot.
    scale : str {'uniform', 'level'}, optional
        Scale the coefficients using the same scale or independently by
        level.
    ax : Axes, optional
        Matplotlib Axes instance

    Returns
    -------
    Figure : Matplotlib figure instance
        Either the parent figure of `ax` or a new pyplot.Figure instance if
        `ax` is None.
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, axisbg='lightgrey')
    else:
        fig = ax.figure

    n_levels = len(coefs)
    n = 2**(n_levels - 1)  # assumes periodic

    if scale == 'uniform':
        biggest = [np.max(np.abs(np.hstack(coefs)))] * n_levels
    else:
        # multiply by 2 so the highest bars only take up .5
        biggest = [np.max(np.abs(i))*2 for i in coefs]

    for i in range(first, n_levels):
        x = np.linspace(2**(n_levels - 2 - i), n - 2**(n_levels - 2 - i), 2**i)
        ymin = n_levels - i - 1 + first
        yheight = coefs[i]/biggest[i]
        ymax = yheight + ymin
        ax.vlines(x, ymin, ymax, linewidth=1.1)

    ax.set_xlim(0, n)
    ax.set_ylim(first - 1, n_levels)
    ax.yaxis.set_ticks(np.arange(n_levels-1, first-1, -1))
    ax.yaxis.set_ticklabels(np.arange(first, n_levels))
    ax.tick_params(top=False, right=False, direction='out', pad=6)
    ax.set_ylabel("Levels", fontsize=14)
    ax.grid(True, alpha=.85, color='white', axis='y', linestyle='-')
    ax.set_title('Wavelet Detail Coefficients', fontsize=16,
                 position=(.5, 1.05))
    fig.subplots_adjust(top=.89)

    return fig


def visualize(data, wavelet, sigma=None):
    coefs = pywt.wavedec(data, wavelet)

    if sigma is None:
        sigma = stand_mad(coefs[-1])

    thresh = sigma*np.sqrt(2*np.log(len(data)))

    denoised = coefs[:]
    denoised[1:] = (pywt.thresholding.soft(i, value=thresh) for i in
                    denoised[1:])

    rec = pywt.waverec(denoised, wavelet)

    plt.plot(data, 'r')
    plt.plot(rec, 'g')
    plt.show()

    return rec


def pdf(data, scale):

    points = list()
    pdf = list()
    dens = list()

    for i in xrange(data.shape[0]):
        dens.append(sm.nonparametric.KDEMultivariate(data[i].T,
                                                     var_type='c',
                                                     bw='cv_ml'))
        points.append(np.linspace(data[i].min(), data[i].max(), 1000))
        pdf.append(dens[i].pdf(points[i]))
        plt.subplot(2, 2, i+1)
        plt.plot(points[i], pdf[i])

    plt.show()

    jdens = sm.nonparametric.KDEMultivariate(data=np.delete(data, [1, 2], 0).T,
                                             var_type='cc',
                                             bw='cv_ml')

    mi = []

    for i in xrange(0, data.shape[1], scale):
        cur_mi = 0
        sample = data[:, i:min(i+scale, data.shape[1])]
        sample = np.delete(sample, [1, 2], 0)
        min_val = np.min(sample)
        max_val = np.max(sample)

        points = np.linspace(min_val, max_val, 1000)
        pt = np.zeros((1, 2))
        for x in points:
            for y in points:
                pt[0][0] = x
                pt[0][1] = y
                jp = jdens.pdf(pt)
                cur_mi += jp*math.log(jp/(dens[0].pdf(x)*dens[2].pdf(y)))

        mi.append(cur_mi)

    plt.plot(mi)
