import numpy as np
import math
from scipy.optimize import minimize
import parse
import matplotlib.pyplot as plt

class Kalman:
    """An implementation of the linear Kalman Filter"""
    dim = None #no. of state variables
    x = None #state
    F = None #state transition function
    G = None #measurement function
    E = None #state-covariance
    R = None #process noise
    Q = None #measurement noise
    x_post = None #posterior states
    E_post = None #posterior covariance

    def __init__(self,R=0.03354,Q=0.04):
        self.dim = 1
        self.x = np.zeros(self.dim)
        
        self.R = np.array(R)
        self.Q = np.array(Q)

        self.E = np.array([0.2])

        self.F = np.identity(self.dim)
        self.G = np.identity(self.dim)

        self.x_post = []
        self.E_post = []
        
    def propagateState(self):
        self.x = np.dot(self.F,self.x);
        self.E = np.dot(self.F,np.dot(self.E,(self.F.T))) + self.R
        
    def makeObservation(self,v):
        eps = v - np.dot(self.G,self.x) 
        S = np.dot(self.G,np.dot(self.E,self.G.T)) + self.Q
        if(S.shape[0] == 1):
            K = np.dot(self.E,np.dot(self.G.T,1/S))
        else:
            K = np.dot(self.E,np.dot(self.G.T,np.linalg.inv(S))) 
        
        self.x = self.x + np.dot(K,eps)
        self.x_post.append(self.x)
        self.E = np.dot((np.identity(self.dim) - np.dot(K,self.G)),self.E) 
        self.E_post.append(self.E)

    def getState(self):
        return self.x

    def getStateCovariance(self):
        return self.E

    def getPosteriorStates(self):
        return self.x_post

    def getPosteriorCovs(self):
        return self.E_post
            
class findBest:
    """Find the best parameters for R and Q by minimizing some objective function"""
    obs = None
    Ri = None
    Qi = None

    def __init__(self,obs,Ri,Qi):
        self.obs = obs
        self.Ri = Ri
        self.Qi = Qi
        
    def predLikelihood(self,v):
        print v
        kf = Kalman(v[0],v[1])
        sum = 0 

        for i in range(len(self.obs)):
            kf.propagateState()
            kf.makeObservation(self.obs[i])
            sum = sum + (np.log(2*np.pi*kf.E) + kf.x*math.pow(kf.E,-1)*kf.x)

        return sum

    def predError(self,v):
        print v
        kf = Kalman(v[0],v[1])
        sum = 0

        for i in range(len(self.obs)):
            kf.propagateState()
            kf.makeObservation(self.obs[i])
            sum = sum + math.pow((-kf.x),2)

        return sum
            
    def getBest(self):
        v0 = np.array([self.Ri,self.Qi])
        cons = ({'type':'ineq','fun':lambda x: x[0] - 1e-4},
                {'type':'ineq','fun':lambda x: x[1] - 1e-4})
        bnds = ((1e-5,1),(1e-5,1))
        res = minimize(self.predLikelihood,v0,constraints=cons,method='COBYLA',options={'disp':True})
        print "Solution is:\n",res.x

if __name__=="__main__":
    files = ['big0.csv','big1.csv','big2.csv','big3.csv']

    [obs,t] = parse.separate(files,0,-1)
    obs_ub = np.copy(obs)

    c = ['b','g','r','c','m','y']

    var = '''
    for i in range(obs.shape[0]):
        print "\n\nFor Gyroscope",(i+1),"\n"
        obs_ub[i,:] = obs[i,:] - np.mean(obs[i,:])
        best = findBest(obs_ub[i,0:50000],0.02,0.02)
        best.getBest()'''


    for i in range(obs.shape[0]):
        obs_ub[i,:] = obs[i,:] - np.mean(obs[i,:])

    kf = Kalman(0.03354,0.04)
    for j in range(obs_ub[i,:].shape[0]):
        kf.propagateState()
        kf.makeObservation(np.mean(obs_ub[:,j]))

    plt.figure(1)
    postStates = np.array(kf.getPosteriorStates())
    plt.subplot2grid((2,2),((int)(i/2),i%2))
    plt.plot(postStates,color=c[i])

    plt.figure(2)
    postCovs = np.array(kf.getPosteriorCovs())
    plt.subplot2grid((2,2),((int)(i/2),i%2))
    plt.plot(postCovs,color=c[i])

    plt.show()
        