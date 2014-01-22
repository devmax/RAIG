import numpy as np
import matplotlib.pyplot as plt
import parse
import scipy.stats as stats

windowSize = 40


def centeredMean(data):

    start = iter(data)
    end = iter(data)

    sum = 0.
    count = 0
    idx = 0
    N = len(data)
    for i in xrange((windowSize)+1):
        sum = sum + next(end)
        count = count + 1

    yield (float)(sum/count)

    for idx in xrange(1, N):
        if idx > (windowSize):
            sum = sum - next(start)
            if idx > (N-(windowSize)-1):
                count = count - 1
            else:
                sum = sum + next(end)
        else:
            count = count + 1
            sum = sum + next(end)

        yield (float)(sum/count)


def removeBias(data):
    if windowSize > data.shape[1]:
        raise Exception("Window size too big!")

    dataf = np.zeros_like(data)
    biasf = np.zeros_like(data)
    #datap = np.zeros_like(data)
    #biasp = np.zeros_like(data)

    for j in xrange(data.shape[0]):
        n = data.shape[1]
        row = data[j, :]
        for i, bias in zip(range(n), centeredMean(row)):
            biasf[j, i] = bias
            dataf[j, i] = data[j, i] - biasf[j, i]

            #biasp[j,i] = np.mean(data[j,max(0,i-windowSize):i+1])
            #datap[j,i] = data[j,i] - biasp[j,i]

    return biasf, dataf  # ,biasp,datap
if __name__ == "__main__":

    files = ['big0.csv', 'big1.csv', 'big2.csv', 'big3.csv']
    [observations, temp] = parse.separate(files)
    for obs in observations:
        [biasf, dataf] = removeBias(obs)

        plt.close('all')

        for i in range(obs.shape[0]):
            biasrate = np.array(np.diff(biasf[i, :]))

            print "\nTest for biasrate is ", stats.kurtosistest(biasrate)
            print "Test for white noise is ", stats.kurtosistest(dataf[i, :])
            print "Test for bias is ", stats.kurtosistest(biasf[i, :])
            print "Test for observation is ",
            stats.kurtosistest(obs[i, :]), "\n"

            plt.figure(i + 1)
            plt.clf()

            plt.subplot2grid((2, 2), (0, 0))
            plt.hist(obs[i, :], color='r')

            plt.subplot2grid((2, 2), (0, 1))
            plt.hist(dataf[i, :], color='g')

            plt.subplot2grid((2, 2), (1, 0))
            plt.plot(dataf[i, :], color='g', label='rate')
            plt.plot(biasf[i, :], color='b', label='bias')
            plt.title('Taking future observations')

            plt.subplot2grid((2, 2), (1, 1))
            plt.hist(biasrate, color='b')

            var = '''
            plt.subplot2grid((2,2),(1,1))
            plt.plot(datap,color='g',label='rate')
            plt.plot(biasp,color='b',label='bias')
            plt.title('Taking only past observations') '''

        plt.figure(5)
        plt.scatter(obs[0, :], obs[1, :])

        plt.figure(6)
        plt.scatter(obs[1, :], obs[2, :])

        plt.figure(7)
        plt.scatter(obs[2, :], obs[3, :])

        plt.figure(8)
        plt.scatter(obs[0, :], obs[3, :])

        plt.show()

        print '*'*150
        print "Covariance in rate is\n", np.corrcoef(dataf)
        print '*'*50
        print "Covariance in bias is\n", np.corrcoef(biasf)
        print '*'*150

        s = raw_input("Press enter to move to next axis!")
