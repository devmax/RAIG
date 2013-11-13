import cv2
import numpy as np
import matplotlib.pyplot as plt
import csv
import parse

def runEM(files=['big0.csv','big1.csv','big2.csv','big3.csv'],scale=0.0001309,skip=0,samples=25000,numClusters=1,covtype=cv2.EM_COV_MAT_GENERIC):

    files = ['big0.csv','big1.csv','big2.csv','big3.csv']

    observations = parse.separate(files,skip,samples)

    print "Observations is a",observations.shape,"type matrix"
        
    em = cv2.EM(numClusters,covMatType=cv2.EM_COV_MAT_GENERIC)
    
    result = em.train(observations.T)
    covMat = em.getMatVector("covs")
    means = em.getMat("means")
    
    print "There are",em.getInt('nclusters'),'clusters'
    print "Means are",means
    print "Weights are",em.getMat("weights")
    print "Covariance Matrix is",covMat

    return means,covMat
    #for arr in covMat:
    #print "Standard deviations for distribution are", arr[0][0]**0.5,arr[1][1]**0.5,arr[2][2]**0.5

if __name__ == "__main__":
    runEM();