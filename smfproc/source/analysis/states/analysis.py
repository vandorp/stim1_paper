""" State analyses """


# System packages
import numpy as np

# Program modules
from smfproc.source.io.datatypes import TimeSeries

# Program modules - local relative imports
from .hmm.hmm import get_hmm, get_best_hmm


### NB: Functions in this module should no longer use the TimeSeries class


#def get_weighted_std(trace, states, descr=None):
#    data = get_state_data(trace, states, descr)
#    data = np.sum([d.data for d in data[0]], axis=0)
#    model = get_state_model(trace, states, descr)
#    return np.std(data-model.data)


def get_snr(data, states, codes, minCnt=3, func=np.min):
    sttMu = get_state_specs(data, states, np.mean, codes)
    sttStd = get_state_specs(data, states, np.std, codes)
    sttCnt = get_state_specs(data, states, len, codes)
    if len(sttMu) == 0:
        return np.nan
    elif len(sttMu) == 1:
        if sttCnt[0] >= minCnt:
            return sttMu/sttStd
        else:
            return np.nan
    else:
        sttSnr = sttMu[sttCnt>=minCnt]/sttStd[sttCnt>=minCnt]
        if len(sttSnr) > 0:
            return func(sttSnr)
        else:
            return np.nan


def get_gamma(don, acc, states, stateCodes):
    donMu = get_state_specs(don, states, np.mean, stateCodes)
    accMu = get_state_specs(acc, states, np.mean, stateCodes)

    donDiff = donMu[1]-donMu[0]
    accDiff = accMu[0]-accMu[1]
    if not (donDiff==0):
        gamma = accDiff/donDiff
    else:
        gamma = np.inf
    return gamma


def get_state_sequence(states):
    if len(states) == 0:
        return np.array([]), np.array([])

    states = np.reshape(states, states.size)
    transitions = np.abs(np.hstack(([1], np.diff(states))))>0
    index = np.arange(states.size)
    return states[transitions], index[transitions]
    

def get_state_data(trace, states, descr):
    # NB this works only for TimeSeries with single trace
    if not (trace.data.shape == states.data.shape):
        print('Error. State data does not match trace data.')
        return []
    
    nrTraces = trace.get_nr_traces()
    data = []
    for n in range(nrTraces):
        stt = states.get(n)
        trc = trace.get(n)
        
        crits = stt.get_crit(descr=descr)
        
        data.append([trc.slice(crit) for crit in crits])
    
    #for c in codes:
    #    crit = states.data == c
    #    data.append(trace.data[:,crit])
    
    #trc = trace.data
    #stt = state.data
    #data = [trc[:, stt==c] for c in codes]

    return data


def get_state_specs(trace, states, func, codes, weighted=False):
    nStates = len(codes)

    data = [trace[states==c] for c in codes]
    
    cnt = np.array([d.size for d in data])
    
    if weighted:
        specs = []
        for n in range(nStates):
            spec = func(data[n])*cnt[n]
            specs.append(spec)
        return np.array(specs)/np.sum(cnt)
    else:
        return np.array([func(d) for d in data])


def get_state_model(trace, states, descr=None):
    # NB this works only for TimeSeries with single trace
    data = get_state_data(trace, states, descr)
    avg = [np.mean(d.data) for d in data[0]]
    crit = states.get_crit(descr=descr)

    if len(crit)==0:
        return TimeSeries([])
    else:
        models = [avg[n]*crit[n] for n in range(len(avg))]
        model = np.sum(np.array(models), axis=0)
        crit = np.sum(np.array(crit), axis=0)==1
        
        crit = np.reshape(crit, crit.size)
        model = np.reshape(model, model.size)
        
        return TimeSeries(model[crit], states.time[crit])


def get_state_transitions(states, descr):
    # NB this works only for TimeSeries with single trace
    if type(descr)==list:
        crit1 = states.get_crit(descr=descr[0])
        crit1 = np.sum(np.array(crit1), axis=0)
        crit2 = states.get_crit(descr=descr[1])
        crit2 = np.sum(np.array(crit2), axis=0)
        
        stt = crit1 + 3*crit2
        crit = np.abs(np.diff(stt)) == 2
    else:
        # NB: This should return two-way transitions between all separate states
        crit = states.get_crit(descr=descr)
        crit = np.squeeze(np.sum(np.vstack(crit), axis=0))
        
        stt = states.data*crit
        crit = np.abs(np.diff(stt)) > 0
        crit = np.squeeze(crit)
    return np.hstack((0,crit))


def get_trace_hmm(data, states, descr, setNrStates, fixNrStates=False):
    # NB this works only for TimeSeries with single trace

    #donSignal = get_state_data(don, states, descr='signal')
    #accSignal = get_state_data(acc, states, descr='signal')

    #data = np.transpose(np.vstack((accSignal[0][0].data, donSignal[0][0].data)))
    stateCrit = states.get_all_crit(descr)
    stateData = []
    for n in range(len(data)):
        #A = get_state_data(data[n], states, descr=descr)
        #stateCrit = states.get_crit(descr=descr)
        #stateCrit = np.sum(np.array(stateCrit), axis=0)
        if np.sum(stateCrit) > 0:
            #stData = data[n].data[0,stateCrit]
            stData = data[n].slice(stateCrit)
            stateData.append(stData.data)
        else:
            return states, None
    
    #stateData = [get_state_data(d, states, descr=descr)[0][0].data for d in data]
    if len(stateData) > 1:
        stateData = np.vstack(stateData)
    else:
        stateData = stateData[0]

    stateData = np.transpose(stateData)
    if fixNrStates==False:
        params = get_best_hmm(stateData, setNrStates)
    else:
        params, success = get_hmm(stateData, setNrStates)
    
    if params == None:
        return states, params

    P = params.postFit
    
    nStates = P.shape[0]
    nPoints = P.shape[1]

    # Collapse probabilities
    stateIndex = np.argmax(P, 0)
    S = np.zeros((nStates, nPoints))
    for n in range(nPoints):
        S[stateIndex[n], n] = 1
    
    # Make state codes
    #codeBase = states.codes['signal']
    #for n in range(nStates):
    #    S[n,:] *= codeBase+n
    #S = np.sum(S, 0)
    #crit = []
    #for n in range(nStates):
    #    crit.append(S[n,:])
    
    # Replace original state with refined states
    #origCrit = states.get_crit(descr=descr)
    origCrit = stateCrit
    #origCrit = np.sum(np.array(origCrit), axis=0)
    origCrit = np.reshape(origCrit, origCrit.size)
    newCrit = np.zeros((nStates, origCrit.size))

    for n in range(nStates):
        newCrit[n,origCrit] = S[n,:]

    crit = [newCrit[n,:] for n in range(nStates)]
    states.set_crit(0, crit, descr)

    return states, params


