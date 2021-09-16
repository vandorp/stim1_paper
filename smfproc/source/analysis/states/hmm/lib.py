import numpy as np
from scipy.cluster.vq import vq, kmeans, whiten
from scipy.stats import norm


length = lambda x:np.max(x.shape)
tp = np.transpose

def TrainPostDec_(p):
    
    #discStates = np.ones(p.nStates)
    x = p.data

    # parameters not loaded from params struct
    relaxiter = 5 # max iter up to which relax A matrix added to trained A matrix

    # initial round of fitting
    A = GetAutoAinitial_(p.nStates, length(x))
    E = GetAutoEinitial_(x, p.nStates, p.nChannels)

    #E[:,:,0] = np.array([[4859,1534],[11265,27645],[35760,41056]])
    #E[:,:,1] = np.array([[2296,3738],[2620,4122],[5386,5543]])

    E = np.array(E, dtype=float)

    A, E, _ = SortStates_(A, E)

    P, F, B, normF, normB, logPx = Post_(A, E, x)
   
    cumsumlognormF = np.cumsum(np.log(normF))
    cumsumlognormB = np.cumsum(np.log(normB[::-1]))[::-1]

    logPxPrevPrev = logPxPrev = logPx #log(P(x|model))
    logPxPermMax = logPx #store logPx to find maximum over constraints permutations

    ## Expectation Maximization of model parameters
    #discStatesList = discStates
    Aperm = np.array(A)
    Eperm = np.array(E)

    everyPermFitCrashed = True
    timedOutIterationNumber = 0


    try:
        if p.nStates==1:
            done = True
        else:
            done = False #done when converged or maxIterEM number of iterations reached

        iterOn = 1
        while (done == False):
            Astore = np.array(A)
            Estore = np.array(E)
            logPxStore = logPx
            Pstore = np.array(P)

            iterOn += 1

            # === reestimate transition rates === #
            Ahat = A*0.0
            e = Get_e_(E,x)
            for i in range(p.nStates):
                for j in range(p.nStates):
                    Ahat[i,j] = \
                            np.sum( \
                            F[i,:-1]*B[j,1:]*A[i,j]*tp(e[1:,j])* \
                            np.exp(cumsumlognormF[:-1]+cumsumlognormB[1:] \
                            -cumsumlognormF[-1]))
            A = np.array(Ahat)

            # here add to A rates up to some iteration 
            if iterOn < relaxiter:
                A = A + (1.0*np.ones((p.nStates,p.nStates))/(1.0*p.nStates))/((1.0*iterOn)**2)
            
            # renormalize A
            A = A/np.reshape(np.repeat(np.sum(A,1), p.nStates), A.shape)
            
            # === reestimate emission parameters === #
            Ep = np.array(E)
            for i in range(p.nStates):
                for j in range(p.nChannels):
                    Ep[i,j,0] = np.sum(P[i,:]*tp(x[:,j]))/np.sum(P[i,:])
                    Ep[i,j,1] = \
                            np.sqrt(np.sum(P[i,:]*((tp(x[:,j]) - Ep[i,j,0])**2)) \
                            /np.sum(P[i,:]))
            E = np.array(Ep)

            # === sort states by mean of emissions in channel 1 ===
            [A,E,rearrangingHappened] = SortStates_(A, E)

            if rearrangingHappened:
                #unzero entries of A that are too small since now different rates in noHops are zeroed
                r_zero, c_zero = np.where(A == 0)
                for r_zero_ix in range(length(r_zero)):
                    A[r_zero[r_zero_ix],c_zero[r_zero_ix]] = 2.0/length(x)
            
            
                       
            try:
                P, F, B, normF, normB, logPx = Post_(A, E, x)
            except:
                #=== HACK! ===
                E = E+0.01
                P, F, B, normF, normB, logPx = Post_(A, E, x)
                #=============



            cumsumlognormF = np.cumsum(np.log(normF))
            cumsumlognormB = np.cumsum(np.log(normB[::-1]))[::-1]

            #print iterOn, logPxPrev, logPx
            #print iterOn, logPx
            #print logPx, logPxPrev, logPx-logPxPrev
            #print np.abs(np.exp(logPx - logPxPrev) - 1)
            
            if (logPx - logPxPrev) < 100:
                crit1 = np.abs(np.exp(logPx - logPxPrev) - 1) < p.threshEMToConverge
            else:
                crit1 = False
            
            #crit2 = np.abs(np.exp(logPx - logPxPrevPrev) - 1) < p.threshEMToConverge
            
            # === Check if not converging === #
            #if crit2 & (iterOn >= relaxiter):
            if logPx - logPxPrev <= 0:
                #print logPx, logPxPrev, logPx-logPxPrev
                done = True
            else:
                #logPxPrevPrev = logPxPrev
                logPxPrev = logPx

            # === see if should stop iterating because converged === #
            if crit1  & (iterOn >= relaxiter):
                done = True
            else:
                logPxPrev = logPx
            
            if (done == False) & (iterOn >= p.maxIterEM) & (iterOn >= relaxiter):
                done = True

        everyPermFitCrashed = False 
    
    except:
        #print ME
        print('boe!')
        return False

    # === see if improved over previous permutation of discStates and noHops === #
    """
    if logPxPermMax < logPx:
        logPxPermMax = logPx
        #bestDiscStates = discStates;
        
        Aperm = A
        Eperm = E
        
        if iter == p.maxIterEM:
            timedOutIterationNumber = True
        else:
            timedOutIterationNumber = False

    logPxMax = logPxPermMax  #maximum log likelihood of data found
    A = Aperm
    E = Eperm
    """

    if p.nStates > 1:
        A = Astore
        E = Estore
        logPx = logPxStore
        P = Pstore

    # === Get BIC === #
    BIC, freeParams = GetBIC_(logPx, length(x), p.nStates, p.nChannels)
    SNRMat = GetSNRMat_(E)
    
    # === return fit === #
    #P = Post_(A,E,x)
    if p.nStates==1:
        P = np.reshape(P, [1,P.size])
    
    p.A = A
    p.E = E
    p.BIC = BIC
    p.SNRMat = SNRMat
    p.postFit = P
    
    return True


def GetAutoAinitial_(nStates, N):
    # space rates order of magnitude apart
    # min rate set by length of data
    if nStates == 1:
        A = np.array(1)
        A = np.reshape(A, A.size)
        return A

    rmin = 10.0/N
    rmax = 0.1
    rates = np.logspace(np.log10(rmin), np.log10(rmax), nStates)
    rates = np.reshape(rates, [nStates, 1])
   
    A = np.repeat(rates, nStates, axis=1)/(nStates-1)
    
    rates = np.reshape(rates, rates.size)
    A -= np.diag(np.diag(A)) - np.diag(1-rates)

    return A


def kmeans_(x, nClusters, nIter=10):
    # Center zero, divide by std
    w = whiten(x)

    # Cluster
    centroids,_ = kmeans(w, nClusters, nIter)
    C = np.sort(centroids)
    clusters,_ = vq(w, C)

    return clusters


def GetAutoEinitial_(x, nStates, nChannels):
    # Gaussian noise model for each channel
    dataMeans = np.mean(x, axis=0)
    dataStdevs = np.std(x, axis=0)
    
    # Initially each state gets mean and std of the full channel
    E = np.zeros((nStates, nChannels, 2))
    for i in range(nStates):
        for j in range(nChannels):
            m = dataMeans[j]
            if nStates > 1:
                m += dataStdevs[j]*(2*i/(nStates-1) - 1)
            v = dataStdevs[j]
            E[i, j, :] = [m, v]

    # Enforce some discrete states have identically distributed emissions
    #if discStates.size > 0:
    #    E = EnforceDiscStates_(E, discStates, np.ones((nStates, x.shape[0])))

    # Try kmeans clustering
    if nStates > 1:
        try:
            for j in range(nChannels):
                idx = kmeans_(x[:,j], nStates)
                for i in range(nStates):
                    clus = x[idx==i, j]
                    E[i, j, :] = [np.mean(clus), np.std(clus)]
        except:
            print('Unable to find initial emissions parameters guess')
    
    return E

"""
def EnforceDiscStates_(E, discStates, P):
    # Enforce some discrete states have identically distributed emissions

    nChannels = E.shape[1]

    tempState = 0
    for i in range(discStates.size):
        for j in range(nChannels):
            temp = E[tempState, j, :] * 0.0

            for k in range(int(discStates[i])):
                Etemp =  E[tempState+k, j, :]*np.sum(P[tempState+k, :])
                temp += Etemp
            
            for k in range(int(discStates[i])):
                Etemp = temp/np.sum(np.sum(P[tempState:tempState+discStates[i], :]))
                E[tempState+k, j, :] = Etemp
        
        tempState += int(discStates[i])
    
    return E
"""

def SortStates_(A, E):
    # Reorders transition matrix A and emissions matrix E to sort states by
    # increasing channel 1 mean (param 1)

    nStates = E.shape[0]
    #if nStates==1:
    #    return A, E, False

    nChannels = E.shape[1]
    means = E[:, 0, 0]

    ix = np.argsort(means)
    means = means[ix]
    
    if np.min(ix == np.arange(nStates)) == 0: #means currently not ascending
        rearrangingHappened = True

        # generate permutation matrix from ix
        p = np.eye(nStates)
        p = p[ix,:]

        # permute transition matrix
        A = np.inner(np.inner(p, A), np.transpose(p))

        # permute emissions matrix
        Ep = np.array(E)
        for i in range(nStates):
            for j in range(nChannels):
                Ep[i, j, :] = E[ix[i], j, :]
        E = np.array(Ep)
    else:
        rearrangingHappened = False

    return A, E, rearrangingHappened


def Post_(A, E, x):
    # computes forward probability matrix F for HMM
    # with transition matrix A, mean vector m,
    # variance vector v, initial probabilities mu0
    # string x

    # compute mu0 from stationary distribution of A
    if (A.size > 1):
        [D, v] = np.linalg.eig(tp(A))
        v = np.abs(v[:,0])
        D = D[0]
        mu0 = v/np.sum(v)
    else:
        mu0 = 1.0
    
    # initialize F
    F = np.zeros((A.shape[0], length(x)))
    normF = np.zeros(length(x))

    # initialize B
    B = np.zeros((A.shape[0], length(x)))
    normB = np.zeros(length(x))

    e = Get_e_(E, x)

    F[:,0] = mu0*e[0,:]
    normF[0] = np.sum(F[:,0])
    #if normF[0] < 0:
    #    normF[0] = 1.0
    #print normF[0]

    F[:,0] = F[:,0]/normF[0] #normalize to avoid underflow
    
    B[:,-1] = np.ones(A.shape[0])
    normB[-1] = np.sum(B[:,-1])
    B[:,-1] = B[:,-1]/normB[-1] #normalize to avoid underflow

    for i in range(1,length(x)):
        #bi = length(x)-i+1
        F[:,i] = np.dot(tp(A), F[:,i-1]) * tp(e[i,:])
        B[:,-i-1] = np.dot(A, B[:,-i] * tp(e[-i,:]))

        #normalize to avoid underflow
        normF[i] = np.sum(F[:,i])
        F[:,i] = F[:,i]/normF[i]
        
        normB[-i-1] = np.sum(B[:,-i-1])
        B[:,-i-1] = B[:,-i-1]/normB[-i-1]

    P = F*B

    Psum = np.sum(P, 0)
    Psum = np.reshape(Psum, [1, Psum.size])

    P = P/np.repeat(Psum, A.shape[0], 0)
    normF[np.isnan(normF)] = 1.0
    normF[normF<0] = 1.0
    #print normF[0]
    #if normF[0] < 0:
    #    normF[0] = 1.0

    logPx = np.sum(np.log(normF))

    #print np.sum(normF>0)
    #print np.where(normF<0)

    return P, F, B, normF, normB, logPx


def Get_e_(E, x):
    nStates = E.shape[0]
    nChannels = E.shape[1]
    
    M = np.zeros((length(x), nChannels, nStates))
    S = np.zeros((length(x), nChannels, nStates))
    xr = np.zeros((x.shape[0], nChannels, nStates))
    
    """
    for i in range(nStates):
        M[:,0,i] = np.ones(length(x))*E[i,0,0]
        M[:,1,i] = np.ones(length(x))*E[i,1,0]

        S[:,0,i] = np.ones(length(x))*E[i,0,1]
        S[:,1,i] = np.ones(length(x))*E[i,1,1]

        xr[:,:,i] = x
    """

    for i in range(nStates):
        for j in range(nChannels):
            M[:,j,i] = np.ones(length(x))*E[i,j,0]
            S[:,j,i] = np.ones(length(x))*E[i,j,1]
        xr[:,:,i] = x
    
    e = np.prod(norm.pdf(xr, M, S), 1)

    # check to see if numerical underflow occurred
    if np.max(np.max(np.isnan(e))) == 1:
        raise ValueError('numerical underflow occurred when computing emission probabilities')

    if np.min(np.sum(np.abs(e),1)) == 0:
        vzeros = np.sum(np.abs(e),1)==0
        if np.sum(vzeros) < x.shape[0]/10.0:
            r = np.where(vzeros == 1)
            for i in range(length(r)):
                e[r[i],:] = np.ones((nStates,1))/nStates
        
            print('warning: numerical underflow occured when computing emission probabilities')
            print('at %d of %d observations'%(np.sum(vzeros), x.shape[0]))
            print('assuming uniform distribution at those points...')
        else:
            raise ValueError('numerical underflow occurred when computing emission probabilities')

    return e


def GetBIC_(logPxMax, N, nStates, nChannels):
    freeParams = nStates**2 - nStates + 2*nStates*nChannels
    BIC = -2*logPxMax + freeParams * np.log(N)
    return BIC, freeParams


def GetSNRMat_(E):
    nStates = E.shape[0]
    nChannels = E.shape[1]

    SNRMat = np.zeros((nStates,nStates))

    meanMat = np.zeros((nStates,nChannels))
    stdevMat = np.zeros((nStates,nChannels))

    for k in range(nStates):
        for i in range(nChannels):
            meanMat[k,i] = E[k,i,0]
            stdevMat[k,i] = E[k,i,1]

    for k1 in range(nStates-1):
        for k2 in range(k1+1,nStates):
            deltaMean = meanMat[k1,:]-meanMat[k2,:]
            SNRMat[k1,k2] = 0.5*(np.sqrt(np.sum((deltaMean/stdevMat[k1,:])**2)) + \
                    np.sqrt(np.sum((deltaMean/stdevMat[k2,:])**2)))
            SNRMat[k2,k1] = SNRMat[k1,k2]

    return SNRMat
