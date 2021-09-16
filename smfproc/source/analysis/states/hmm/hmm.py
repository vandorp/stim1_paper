""" Hidden Markov Model fitting.

HMM fitting routines were taken from SMART package (Greenfeld, Herschlag PLOS One 2012)
"""

import numpy as np

from .lib import TrainPostDec_

class Container():
    pass


def get_hmm(data, nStates):
    """ Fit HMM """

    params = Container()

    # Maximum number of iterations for EM to converge
    params.maxIterEM = 200

    # likelihood threshold for EM to converge, don't make lower than this to
    # avoid strange behavior
    params.threshEMToConverge = 10**(-4)

    # Warns if 2 states have SNR less than this
    params.SNRwarnthresh = 1
    
    # Number of channels - acceptor is channel 0, donor is channel 1
    params.nChannels = data.shape[1]

    # Number of states in the model
    params.nStates = nStates

    # x is T by num_channels, column 1 is acceptor (red), column 2 is donor 
    # (green)
    params.data = data;
    
    success = TrainPostDec_(params)

    return params, success

"""
def mu_diff_too_small(mu):
    print np.diff(mu)
    print np.mean(mu)
    return np.sum(np.diff(mu) < MIN_STATE_DIFF*np.mean(mu)) > 0
"""
def get_best_hmm(data, maxNrStates):
    """ Implement BIC to get best model. """
    
    try:
        p, success = get_hmm(data, 1)
    except:
        #print 'Error at state 1'
        return None

    for n in range(1, maxNrStates):

        #print n, p.BIC
        try:
            pNew, success = get_hmm(data, n+1)
            if not success:
                pNew = p
                continue
            #return p
        except:
            #print 'Error at state %d'%(n+1)
            pNew = p
            continue
        
        #print pNew.BIC, p.BIC
        if pNew.BIC >= 1.00*p.BIC:
            #print pNew.BIC
            #print 'Minimum at state %d'%(n+1)
            return p #, n
            #p = pNew
        else:
            p = pNew
    
    #print 'Stopped at max nr states'
    return pNew #, n+1


"""
def get_best_hmm(data):
    
    allP = []

    try:
        p = get_hmm(data, 1)
        allP.append(p)
    except:
        return None

    for n in range(1, maxNrStates):

        try:
            p = get_hmm(data, n+1)
            allP.append(p)
        except:
            break
            #return allP[n-1]

        #print pNew.BIC, p.BIC
        if allP[n].BIC >= allP[n-1].BIC:
            del allP[-1]
            break
            #return allP[n-1] #, n
        #else:
        #    p = pNew

    statesTooClose = [False]
    for m in range(1, n):
        #accStatesTooClose = mu_diff_too_small(allP[m].E[:,0,0])
        #donStatesTooClose = mu_diff_too_small(allP[m].E[:,1,0])
        #statesTooClose = mu_diff_too_small(np.sum(allP[m].E[:,:,0], axis=1))

        #print accStatesTooClose
        #print statesTooClose

        E = allP[m].E

        accDiff = np.diff(E[:,0,0])
        donDiff = np.diff(E[:,1,0])
        muDiff = np.diff(E[:,:,0], axis=0)
        muRef = np.mean(np.sum(E[:,:,0], axis=1))

        crit = np.sum(muDiff > MIN_STATE_DIFF*muRef) == 0
        statesTooClose.append(crit)
    
    index = np.array(statesTooClose) == False
    index = np.max(np.arange(n)[index])

    return allP[index] #, n+1
"""
