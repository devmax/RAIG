import numpy as np
import matplotlib.pyplot as plt
import random
import csv
import cv2
import runEM
import parse

class KalmanFilter:
    """An implementation of the Linear Kalman Filter"""
    dim = None
    h = None
    F = None
    Eh = None
    Ev = None
    A = None
    B = None
    h_post = None
    F_post = None
    F_pre = None
    
    def __init__(self,h0_em=None,Ev_em=None):
        self.dim = 4 #number of states
        self.h = np.zeros(self.dim)
        if(h0_em != None):
            self.h = h0_em
        
        rateProcNoise = 0.03354 #std. dev for proc. noise
        biasProcNoise = rateProcNoise/45 #std. dev of bias noise
        measNoise = 0.04
        initRateNoise = 0.000001
        initBiasNoise = 0.1

        self.F = np.identity(self.dim)
        self.F = self.F * (initBiasNoise*initBiasNoise)
        #self.F[0,0] = initRateNoise*initRateNoise

        self.Eh = np.identity(self.dim)
        self.Eh = self.Eh * biasProcNoise*biasProcNoise
        #self.Eh = self.Eh * (biasProcNoise*biasProcNoise)
        #self.Eh = self.Eh + 0.00001
        #self.Eh[0,1:self.dim] = 0
        #self.Eh[1:self.dim,0] = 0
        #self.Eh[0,0] = rateProcNoise*rateProcNoise
        
        if(Ev_em == None):
            self.Ev = (np.identity(self.dim)) * (measNoise*measNoise)
        else:
            self.Ev = Ev_em
            
        self.A = (np.identity(self.dim))

        self.B = np.array([[1,1,0,0,0],[1,0,1,0,0],[1,0,0,1,0],[1,0,0,0,1]])
        self.B = np.identity(self.dim)

        self.h_post = []
        self.F_post = []
        self.F_pre = []

    def propagateState(self):
        self.h = np.dot(self.A,self.h);
        self.F = np.dot(self.A,np.dot(self.F,(self.A.T))) + self.Eh
        self.F_pre.append(self.F)

    def makeObservation(self,v):
        eps = v - np.dot(self.B,self.h) 
        S = np.dot(self.B,np.dot(self.F,self.B.T)) + self.Ev 
        K = np.dot(self.F,np.dot(self.B.T,np.linalg.inv(S))) 
        
        self.h = self.h + np.dot(K,eps)
        self.h_post.append(self.h)
        self.F = np.dot((np.identity(self.dim) - np.dot(K,self.B)),self.F) 
        self.F_post.append(self.F)
        
    def getState(self):
        return self.h.tolist()

    def getCovariance(self):
        return self.F.tolist()

    def getPosteriorStates(self):
        return self.h_post

    def getPosteriorCovs(self):
        return self.F_post

    def getPriorCovs(self):
        return self.F_pre

class SmoothedKalman:
    hs = None
    Fs = None
    Eh = None
    Ev = None
    A = None
    B = None
    hkf = None
    Fkf = None
    Fkfp = None
    dim = None
    Fsc = None
    
    def __init__(self,kf):
        self.dim = kf.dim
        self.hkf = kf.getPosteriorStates()
        self.Fkfp = kf.getPriorCovs()
        self.Fkf = kf.getPosteriorCovs()
        self.hs = self.hkf[:]
        self.Fs = self.Fkf[:]

        self.Eh = np.copy(kf.Eh)
        self.Ev = np.copy(kf.Ev)
        self.A = np.copy(kf.A)
        self.B = np.copy(kf.B)

        self.Fsc = self.Fkf[:]
        self.Fsc.pop()
        KT =  np.dot(self.Fkfp[-1],np.dot(self.B.T,np.linalg.inv(np.dot(self.B,np.dot(self.Fkfp[-1],self.B.T)) + self.Ev)))
        self.Fsc[-1] = np.dot((np.identity(self.dim) - np.dot(KT,self.B)),np.dot(self.A,self.Fkf[-2]))

        self.Fkfp.append(self.Eh + np.dot(self.A,np.dot(self.Fkf[-1],self.A.T)))

    def backRecurse(self):

        #As = [((self.A*self.Fkf[i]).getT())*((self.A*self.Fkf[i]*(self.A.getT())).getI()) for i in range(len(self.Fkf)) ]
        #S = [(self.Fkf[i] - As[i]*self.A*self.Fkf[i]) for i in range(len(self.Fkf))]
                
        J = [np.dot(self.Fkf[i],np.dot(self.A.T,np.linalg.inv(self.Fkfp[i+1]))) for i in range(len(self.Fkf))]

        for i in range(1,len(self.hkf)):
            self.hs[-(i+1)] = self.hkf[-(i+1)] + np.dot(J[-(i+1)],(self.hs[-i] - np.dot(self.A,self.hkf[-(i+1)])))
            self.Fs[-(i+1)] = self.Fkf[-(i+1)] + np.dot(J[-(i+1)],np.dot((self.Fs[-i] - self.Fkfp[-i]),J[-(i+1)].T)) 

            #self.hs[-(i+1)] = As[-(i+1)]*self.hs[-i] + self.hkf[-(i+1)] - As[-(i+1)]*self.A*self.hkf[-(i+1)]
            #self.Fs[-(i+1)] = As[-(i+1)]*self.Fs[-i]*(self.A.getT()) + S[-(i+1)]
            #self.Fs[-(i+1)] = 0.5*(self.Fs[-(i+1)] + self.Fs[-(i+1)].getT())

            #self.Fsc[-(i+1)] = As[-(i+1)]*self.Fs[-(i+1)] + self.hs[-(i+1)]*(self.hkf[-(i+1)].getT()) 
            if i<(len(self.hkf)-1):
                self.Fsc[-(i+1)] = np.dot(self.Fkf[-(i+1)],J[-(i+2)].T) + np.dot(J[-(i+1)],np.dot((self.Fsc[-i] - np.dot(self.A,self.Fkf[-(i+1)])),J[-(i+2)].T))

    def getSmoothedStates(self):
        return self.hs

    def getSmoothedCovs(self):
        return self.Fs

    def getCrossCovs(self):
        return self.Fsc

class iterateEM:
    T = None
    An = None
    Bn = None
    Ehn = None
    Evn = None
    xon = None
    Fon = None

    hs = None
    Fs = None
    Fsc = None
    v = None
    
    E1 = None
    E2 = None
    E1m1 = None
    Evh = None
    Evv = None
    Ehh = None

    A = None
    B = None

    dim = None
    
    def __init__(self,skf,obs):

        self.dim = skf.dim
        self.hs = skf.hs[:]
        self.Fs = skf.Fs[:]
        self.Fsc = skf.Fsc[:]

        self.v = obs
        self.T = len(self.hs)

        self.A = np.copy(skf.A)
        self.B = np.copy(skf.B)
        
    def calcSums(self):
        self.E1 = [(np.outer(self.hs[i],self.hs[i]) + self.Fs[i]) for i in range(self.T)] #Eznzn
        self.E1 = sum(self.E1) 
        
        self.E2 = [(np.outer(self.hs[i+1],self.hs[i]) + self.Fsc[i]) for i in range(self.T-1)] #Eznznm1
        self.E2 = sum(self.E2)
        
        self.E1mN = self.E1 - (np.outer(self.hs[-1],self.hs[-1]) + self.Fs[-1])
        self.E1m1 = self.E1 - (np.outer(self.hs[0],self.hs[0]) + self.Fs[0])
        
        self.Evh = [np.outer(self.v[:,i],self.hs[i]) for i in range(self.T)] 
        self.Evh = sum(self.Evh)
        
        self.Evv = [np.outer(self.v[:,i],self.v[:,i]) for i in range(self.T)]
        self.Evv = sum(self.Evv)

        self.Ehv = [np.outer(self.hs[i],self.v[:,i]) for i in range(self.T)] 
        self.Ehv = sum(self.Ehv)
        
    def update(self):
        self.calcSums()
        
        #self.An = np.dot(self.E2,np.linalg.inv(self.E1mN))
        self.An = self.A

        self.Ehn = (self.E1m1 - np.dot(self.An,self.E2.T))*(1/(self.T-1))

        #self.Bn = np.dot(self.Evh,np.linalg.inv(self.E1))
        self.Bn = self.B

        self.Evn = (self.Evv - np.dot(self.Bn,self.Evh.T))*(1/self.T)
        self.xon = self.hs[0]
        self.Fon = self.Fs[0]

if __name__ == "__main__":
    numIter = 15
    #trueRate = [1.2,1.7]
    #obsNoise = [0.15,0.15]

    files = ['big0.csv','big1.csv','big2.csv','big3.csv']
    scale = 0.0001309
    samples = 25000
    
    #[h0_em,Ev_em] = runEM.runEM(files,scale,0,samples)

    total = 1020000
    skip = total - (samples+10000)
    samples = 35000
    observations = parse.separate(files,0,samples)
    
    kf = [KalmanFilter()]

    N = observations.shape[1]
    #plt.subplot(2,2,1)
    #plt.hist([o.tolist()[0] for o in observations])
    #plt.subplot(2,2,2)
    #plt.hist([o.tolist()[1] for o in observations])
    #plt.subplot(2,2,3)
    #plt.hist([o.tolist()[2] for o in observations])
    #plt.subplot(2,2,4)
    #plt.hist([o.tolist()[3] for o in observations])

    #plt.show()
    #raw_input("Waiting...")
    
    for i in range(N):
        kf[-1].propagateState()
        kf[-1].makeObservation(observations[:,i])

    skf = [SmoothedKalman(kf[-1])]
    skf[-1].backRecurse()

    posteriorState = kf[-1].getPosteriorStates()
    posteriorCovariance = kf[-1].getPosteriorCovs()

    smoothPostState = skf[-1].getSmoothedStates()
    smoothPostCov = skf[-1].getSmoothedCovs()

    colors = ['b','g','r','m','c']
    
    plt.figure(1)
    plt.clf()
#    plt.subplot(2,3,1)
    plt.subplot(2,1,1)
    plt.xlabel('Iteration number')
    for i in range(posteriorState[0].shape[0]):
        plt.plot([s[i] for s in posteriorState],color=colors[i])

#    plt.subplot(2,3,2)
#    for i in range(observations.shape[0]):
#        plt.plot([observations[i][j] for j in range(observations.shape[1])],color=colors[i])
    
#    plt.subplot(2,3,3)
#    for i in range(posteriorCovariance[0].shape[0]):
#        plt.plot([s[i,i] for s in posteriorCovariance],color=colors[i])

#    plt.subplot(2,3,4)
    plt.subplot(2,1,2)
    plt.xlabel('Iteration number')
    for i in range(smoothPostState[0].shape[0]):
        plt.plot([s[i] for s in smoothPostState],color=colors[i])
    
#    plt.subplot(2,3,5)
#    for i in range(smoothPostCov[0].shape[0]):
#        plt.plot([s[i][i] for s in smoothPostCov],color=colors[i])

    plt.show()

    angles = np.array(posteriorState)[:,0]
    angles = np.cumsum(angles)
    angles = angles / 45

    sangles = np.array(smoothPostState)[:,0]
    sangles = np.cumsum(sangles)
    sangles = (angles / 45)*180/np.pi
    
    plt.figure(2)
    plt.clf()    

    plt.subplot2grid((2,2),(0,0),colspan=2)
    plt.ylabel('Angle in degrees')
    plt.xlabel('Time in minutes')
    plt.xticks(np.arange(0,samples,45*60),np.arange(0,samples/(45*60),1))
    plt.plot(angles,color = colors[0])
    
    plt.subplot2grid((2,2),(1,0))
    plt.ylabel('Observations from various gyros in degrees')
    plt.xlabel('Time in minutes')
    plt.xticks(np.arange(0,samples,45*60),np.arange(0,samples/(45*60),1))
    for i in range(observations.shape[0]):
        plt.plot(observations[0,:],color=colors[i+1])

    plt.subplot2grid((2,2),(1,1))
    plt.ylabel('Average observation in degrees')
    plt.xlabel('Time in minutes')
    avgAngle = np.average(observations,axis=0)
    plt.plot(avgAngle,color=colors[0])

    plt.show()
    
    s = raw_input("Press enter to continue!")
    
    for i in range(numIter):
        em = iterateEM(skf[-1],observations)
        em.update()

        print "Final values after ",i+1," iterations are \n An=",em.An,"\n Eh=",em.Ehn,"\n Bn=",em.Bn,"\n Ev=",em.Evn,"\n x0=",em.xon,"\n F0=",em.Fon

        kf.append(KalmanFilter(em.xon,em.Evn))
        kf[-1].A = em.An
        kf[-1].B = em.Bn
        kf[-1].Eh = em.Ehn
        kf[-1].Ev = em.Evn
        kf[-1].h = em.xon
        kf[-1].F = em.Fon

        for i in range(N):
            kf[-1].propagateState()
            kf[-1].makeObservation(observations[:,i])

        skf.append(SmoothedKalman(kf[-1]))
        skf[-1].backRecurse()

        del em

        plt.clf()
        posteriorState = kf[-1].getPosteriorStates()
        posteriorCovariance = kf[-1].getPosteriorCovs()

        smoothPostState = skf[-1].getSmoothedStates()
        smoothPostCov = skf[-1].getSmoothedCovs()

        plt.figure(1)
        plt.clf()
    #    plt.subplot(2,3,1)
        plt.subplot(2,1,1)
        for i in range(posteriorState[0].shape[0]):
            plt.plot([s[i] for s in posteriorState],color=colors[i])

    #    plt.subplot(2,3,2)
    #    for i in range(observations.shape[0]):
    #        plt.plot([observations[i][j] for j in range(observations.shape[1])],color=colors[i])
            
    #    plt.subplot(2,3,3)
    #    for i in range(posteriorCovariance[0].shape[0]):
    #        plt.plot([s[i,i] for s in posteriorCovariance],color=colors[i])

    #    plt.subplot(2,3,4)
        plt.subplot(2,1,2)
        for i in range(smoothPostState[0].shape[0]):
            plt.plot([s[i] for s in smoothPostState],color=colors[i])
    
    #    plt.subplot(2,3,5)
    #    for i in range(smoothPostCov[0].shape[0]):
    #        plt.plot([s[i][i] for s in smoothPostCov],color=colors[i])
    
        plt.show()

        s = raw_input("Press enter to process to next iteration!")
