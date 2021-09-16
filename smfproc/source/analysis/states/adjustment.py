""" State adjustments"""

# System packages
import numpy as np

# Program modules - local relative imports
from .analysis import get_state_specs


def correct_bleach_time():
    cStates = states.get(n)
    has_acceptor_bleach = cStates.has('acceptor_bleached')
    if has_acceptor_bleach:
        crit = cStates.get_all_crit('signal')
        lastSignalState = cStates.data[crit][-1]

        crit = cStates.get_all_crit('acceptor_bleached')
        firstBleachState = cStates.data[crit][0]

        #donPostMu = get_state_specs(don.get(n).data, cStates.data, 
        #        np.mean, [firstBleachState])
        #donPreMu = get_state_specs(don.get(n).data, cStates.data, 
        #        np.mean, [lastSignalState])
        accPostMu = get_state_specs(acc.get(n).data, cStates.data, 
                np.mean, [firstBleachState])
        accPreMu = get_state_specs(acc.get(n).data, cStates.data, 
                np.mean, [lastSignalState])
        accPostStd = get_state_specs(acc.get(n).data, cStates.data, 
                np.std, [firstBleachState])
        accPreStd = get_state_specs(acc.get(n).data, cStates.data, 
                np.std, [lastSignalState])

        accAvgStd = (accPostStd[0]+accPreStd[0])/2
        
        #crit1 = donPostMu > donPreMu
        crit1 = accPostMu[0] < accPreMu[0]
        crit2 = accPreMu[0]-accPostMu[0] > stateSettings['BLEACH_SNR']*accAvgStd

        if not (crit1 and crit2):
            # acceptor_bleached states will be changed into newly
            # created signal states
            sigCodes = cStates.get_state_codes('signal')
            blCodes = cStates.get_state_codes('acceptor_bleached')
            blWinCode = cStates.codes['acceptor_bleach_event']
            blCodes = blCodes + [blWinCode]

            # if there is one state after acceptor bleach
            newCodes = sigCodes[-1] + np.arange(len(blCodes)) + 1
            for m in range(len(blCodes)):
                # replace codes
                crit = cStates.data==blCodes[m]
                crit = np.reshape(crit, crit.size)
                states.data[n,crit] = newCodes[m]

            # if there are >1 states after acceptor bleach
            # move the acceptor bleach the next state transition. this
            # will have to be a transition to a lower acceptor state
            # because otherwise it would have been detected earlier.



def merge_states(don, acc, states, descr, snr, accountForSign=False):
    for m in range(3):
        statesMerged = merge_similar_states(don, acc, states, descr, snr, accountForSign)
        if statesMerged==False:
            return states.data

"""
            for m in range(3):
                statesMerged = merge_similar_states(
                        don.get(n), acc.get(n), cStates, 
                        'acceptor_bleached', stateSettings['TRANSITION_SNR'])
                
                if statesMerged:
                    states.data[n,:] = cStates.data
                    pass
                else:
                    break
"""



def merge_similar_states(don, acc, states, descr, minDiff, accountForSign):
    
    stt = states.data
    codes = states.get_codes(descr)
    nrStates = len(codes)

    if nrStates < 2:
        return False

    donMu = get_state_specs(don.data, stt, np.mean, codes)
    accMu = get_state_specs(acc.data, stt, np.mean, codes)
    
    # Calculate the average fluorescence as a reference, to decide if
    # a fluctuation is significant
    #cnt = get_state_specs(stt, stt, np.size, codes)
    
    #donAvg = np.sum(np.array(donMu)*np.array(cnt))/np.sum(cnt)
    #accAvg = np.sum(np.array(accMu)*np.array(cnt))/np.sum(cnt)
    #refMu = donAvg + accAvg

    donStd = get_state_specs(don.data, stt, np.std, codes)
    accStd = get_state_specs(acc.data, stt, np.std, codes)

    for m in range(nrStates):
        for k in range(nrStates):
            if not (m==k):
                #donMuDiff = np.abs(donMu[m]-donMu[k])/refMu < minDiff
                #accMuDiff = np.abs(accMu[m]-accMu[k])/refMu < minDiff

                
                donAvgStd = (donStd[m] + donStd[k])/2
                crit1 = np.abs(donMu[m]-donMu[k]) < minDiff*donAvgStd
                accAvgStd = (accStd[m] + accStd[k])/2
                crit2 = np.abs(accMu[m]-accMu[k]) < minDiff*accAvgStd

                donSign = donMu[m]-donMu[k] > 0
                accSign = accMu[m]-accMu[k] > 0

                crit3 = (donSign==accSign)
                
                if accountForSign:
                    crit4 = crit1 and crit2 and crit3
                else:
                    crit4 = crit1 and crit2


                #if (donMuDiff and accMuDiff):
                if crit4:
                    crit = np.logical_or(stt==codes[m], stt==codes[k])
                    
                    # Use the lowest value code
                    newCode = np.min([codes[m], codes[k]])
                    stt[crit] = newCode
                    return True
                    
    return False


"""
def correct_baseline(traces, states):

    # Baseline correction
    nrTraces = traces[0].get_nr_traces()

    donParams = states.get_parameters(traces[0])
    accParams = states.get_parameters(traces[1])

    donMinStates = [p['avg'][p['codes']>=s.STATE_CODES['donBleached']] 
            for p in donParams]
    accMinStates = [p['avg'][p['codes']>=s.STATE_CODES['accBleached']] 
            for p in accParams]
    
    for n in range(nrTraces):
        if len(donMinStates[n]) > 0:
            traces[0].data[n,:] -= donMinStates[n][0]

        if len(accMinStates[n]) > 0:
            traces[1].data[n,:] -= accMinStates[n][0]

def check_zero_states(traces, states):

    nrTraces = traces[0].get_nr_traces()

    donParams = states.get_parameters(traces[0])
    accParams = states.get_parameters(traces[1])

    for n in range(nrTraces):
        a = accParams[n]
        accAvg = a['avg'][a['codes']==s.STATE_CODES['accBleached']]
        accStd = a['std'][a['codes']==s.STATE_CODES['accBleached']] 

        if len(accAvg) > 0:
            accIsZero = accAvg[0] < s.STATE_SNR*accStd[0]
            if not accIsZero:
                inds = states.data[n,:] == s.STATE_CODES['accBleached']
                states.data[n,inds] = s.STATE_CODES['signal']

        d = donParams[n]
        donAvg = d['avg'][d['codes']==s.STATE_CODES['donBleached']]
        donStd = d['std'][d['codes']==s.STATE_CODES['donBleached']] 

        if len(donAvg) > 0:
            donIsZero = donAvg[0] < s.STATE_SNR*donStd[0]
            if not donIsZero:
                inds = states.data[n,:] == s.STATE_CODES['donBleached']
                states.data[n,inds] = s.STATE_CODES['accBleached']

def adjust(traces, states):
    #check_zero_states(traces, states)
    correct_baseline(traces, states)
"""
