from math import sqrt
import random
import matplotlib.pyplot as plt
import numpy as np
from operator import sub 
import cv2,cv

def populateGaussian(trueRate,mu,stdDev,rspec,num,T,dt,s=None):
    rate = []
    for i in range(num):
        if(s != None):
            random.seed(s[i])
        else:
            random.seed()
        rate.append([trueRate + random.gauss(mu[i],stdDev[i]) for x in range(int(T/dt))])
    return rate    

def integrateRate(rate,dt):
    dangle = [[r*dt for r in row] for row in rate]
    angle = [np.cumsum(row) for row in dangle]
    return angle

def averageRate(rate):
    return np.mean(rate,axis=0)

def averageAngle(angle):
    return np.mean(angle,axis=0)

def findLeastError(truth,meas):
    sqErr = []
    
    for row in meas:
        diff = map(sub,truth,row)
        diffSq = [x**2 for x in diff]
        sqErr.append(sum(diffSq))

    return sqErr.index(min(sqErr))

def findBestRate(trueRate,rate):
    trueRate = [0 for i in range(len(rate[0]))]
    return rate[findLeastError(trueRate,rate)]

def findBestAngle(trueRate,angle,dt):
    trueAngle = [trueRate*dt for i in range(len(angle[0]))]
    return angle[findLeastError(trueAngle,angle)]
    
def plotData(numRows,numCols,spc,T,dt,*args):
     for i in range(numRows*numCols):
        plt.subplot(numRows,numCols,i+1)
        plt.xticks(np.arange(0,T/dt,spc/dt),np.arange(0,T/spc,spc/60))
        for row in args[0][i]:
            plt.plot(row)
            
    plt.show()

if __name__ == "__main__":
    dt = 0.02 #sampling time
    T = 180 #total time in seconds
    num = 3 #number of gyros
    stdDev = [0.035,0.02,0.03] #standard deviation of white noise in gyro output
    rspec = 0.05
    trueRate = 0 #ground truth for gyro rate
    s = range(num)
    mu = [0,0,0]
    
    rate = populateGaussian(trueRate,mu,stdDev,rspec,num,T,dt,s)
    angle = integrateRate(rate,dt)

    avgRate = [averageRate(rate)]
    avgAngle = [averageAngle(angle)]

    bestRate = [findBestRate(trueRate,rate)]
    bestAngle = [findBestAngle(trueRate,angle,dt)]

    data = [rate,angle,avgRate,avgAngle,bestRate,bestAngle]

    spc = 60.0

    plt.figure(1)
    plotData(3,2,spc,T,dt,data)
    plt.show()
    
    #em = cv2.EM(covMatType=cv2.EM_COV_MAT_GENERIC)

    #result = em.train(np.transpose(rate))

    #print "There are",em.getInt('nclusters'),'clusters'
    #print "Means are",em.getMat("means")
    #covMat = em.getMatVector("covs")
    #print "Covariance Matrix is",covMat
    #print "Weights are",em.getMat("weights")

    #for arr in covMat:
    #print "Standard deviations for distribution are", arr[0][0]**0.5,arr[1][1]**0.5,arr[2][2]**0.5
    
    #plt.figure(2)
    
    #plt.subplot(3,1,1)
    #plt.fill_between(range(len(linearList)),result[3][:,0],0)

    #plt.subplot(3,1,2)
    #plt.fill_between(range(len(linearList)),result[3][:,1],0)

    #plt.subplot(3,1,3)
    #plt.fill_between(range(len(linearList)),result[3][:,2],0)
    
    #plt.show()

    