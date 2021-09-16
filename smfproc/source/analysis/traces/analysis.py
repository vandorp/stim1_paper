""" Low level trace analysis routines."""

# System packages
import numpy as np

# Program modules - relative imports
from .filters import SavGol, moving_rms, moving_average, MovingAverage

# Program modules
from smfproc.source.io.datatypes import TimeSeries, SmfSelect, State
from smfproc.source.analysis.states.analysis import get_state_transitions



def correct_xtalk(don, acc, xtalk):
    acc.data -= xtalk*don.data


def correct_baseline(don, acc, states):
    nrTraces = don.get_nr_traces()

    hasAccBleach = states.has('acceptor_bleached')
    hasDonBleach = states.has('donor_bleached')

    for n in range(nrTraces):
        if hasAccBleach[n]:
            accBleached = states.get_crit(n, 'acceptor_bleached')[0][0]
            
            baseline = np.mean(acc.data[n, accBleached])
            #noise = np.std(acc.data[n, accBleached])

            #if baseline < ZERO_THRESHOLD*noise:
            acc.data[n,:] = acc.data[n,:] - baseline

        if hasDonBleach[n]:
            donBleached = states.get_crit(n, 'donor_bleached')[0][0]

            baseline = np.mean(acc.data[n, donBleached])
            #noise = np.std(acc.data[n, donBleached])
            #if baseline < ZERO_THRESHOLD*noise:
            acc.data[n,:] = acc.data[n,:] - baseline

            baseline = np.mean(don.data[n, donBleached])
            #noise = np.std(don.data[n, donBleached])
            #if baseline < ZERO_THRESHOLD*noise:
            don.data[n,:] = don.data[n,:] - baseline


def get_fret_ratio(don, acc):
    tot = don.data + acc.data

    # Prevent zero-division
    tot[tot==0] = np.inf
    fret = acc.data/tot
    fret[fret==0] = np.inf

    return TimeSeries(fret, don.time)


def correct_fret(fret, gamma):
    for n in range(len(gamma)):
        G = gamma[n]
        if np.isnan(G):
            fret[n,:] = fret[n,:]*np.nan
        elif G==0:
            fret[n,:] = fret[n,:]*np.nan
        else:
            #denom = fret[n,:]+G-G*fret[n,:]
            fretn = np.array(fret[n,:]*10000, dtype=int)/10000.0
            denom = fretn+G-G*fretn
            fret[n,:] = fret[n,:]/denom
    #gamma = np.reshape(gamma, [gamma.size, 1])
    #G = np.repeat(gamma, fret.shape[1], 1)
    #return fret/(fret+G-G*fret)

def correct_donor(don, gamma):
    for n in range(len(gamma)):
        don[n,:] = don[n,:]*gamma[n]


""" Trace selection """

def get_initial_trace_selection(don, acc, states,
        filterWindow,
        stateSNR):
    
    nrTraces = don.get_nr_traces()
    
    crit = states.data < states.codes['acceptor_present']
    sig = don.data + acc.data
    donSig = sig[:,crit[0,:]]
    rms = np.array([moving_rms(donSig[n,:], filterWindow) 
        for n in range(nrTraces)])
    avg = np.array([moving_average(donSig[n,:], filterWindow) 
        for n in range(nrTraces)])
    #snr = avg/rms


    ####
    
    # Introduce categories:
    #   - acceptor+donor bleach
    #   - acceptor bleach only
    #   - donor bleach only
    #   - no acceptor and no donor bleach
    #   - donor signal only, but acceptor confirmed


    ####

    # --- Criteria --- #
    hasSignal = states.has('signal')
    #hasInitSignal = avg[:,0]/np.median(rms, axis=1) >= stateSNR
    hasAcceptor = states.has('acceptor_present')

    # --- Store selections --- #
    select = SmfSelect()

    crit = np.logical_and(hasSignal==False, hasAcceptor==False)
    #crit = hasSignal==False
    select.set(crit, 'no_fret')
    
    crit = np.logical_and(hasSignal==False, hasAcceptor==True)
    select.set(crit, 'zero_fret')

    #crit = np.logical_and(hasSignal, hasInitSignal==False)
    #select.set(crit, 'no_initial_signal')
    
    return select



""" --- TRACE ANALYSIS ----------------------------------------------------- """

def get_signal_states(allTraces, descr,
        filterWindow,
        savGolOrder,
        traceSnr,
        bleachWindow,
        signalWindow):
   
    laser532 = allTraces[descr.index('laser532')]
    laser638 = allTraces[descr.index('laser638')]
    
    don = allTraces[descr.index('donor')]
    acc = allTraces[descr.index('acceptor')]

    donorLaserOn = laser532.data > 0
    acceptorLaserOn = laser638.data > 0

    donorLaserOn = np.reshape(donorLaserOn, donorLaserOn.size)
    acceptorLaserOn = np.reshape(acceptorLaserOn, acceptorLaserOn.size)
   
    donIllum = [don.slice(donorLaserOn), acc.slice(donorLaserOn)]
    accIllum = [don.slice(acceptorLaserOn), acc.slice(acceptorLaserOn)]    
    donIllumFilt = [donIllum[0].copy(), donIllum[1].copy()]

    smoothFilter = SavGol(filterWindow, savGolOrder)
    smoothFilter.run(donIllumFilt[0])
    smoothFilter.run(donIllumFilt[1])

    nrTraces = don.get_nr_traces()
    nrFrames = don.time.size
    states = State(np.zeros(don.data.shape))

    for n in range(nrTraces):
        
        donRms = np.median(moving_rms(donIllum[0].get(n).data, filterWindow))
        accRms = np.median(moving_rms(donIllum[1].get(n).data, filterWindow))

        donThresh = donRms*traceSnr
        accThresh = accRms*traceSnr

        donZeroThresh = donRms*traceSnr
        accZeroThresh = accRms*traceSnr
        
        #donBleachIndex = get_bleach_index(donIllumFilt[0].get(n).data, donThresh)
        #accBleachIndex = get_bleach_index(donIllumFilt[1].get(n).data, accThresh)

        #print donBleachIndex, accBleachIndex
        #print 'don'
        donBleachIndex = iteratively_get_bleach_index(
                donIllumFilt[0].get(n).data, donThresh, donZeroThresh,
                bleachWindow, signalWindow)
        #print 'acc'
        accBleachIndex = iteratively_get_bleach_index(
                donIllumFilt[1].get(n).data, accThresh, accZeroThresh,
                bleachWindow, signalWindow)
        #print donBleachIndex, accBleachIndex


        states.set_transition(n, 0, 'signal', 0)
        
        if donBleachIndex <= 0:
            # This means there was no donor signal detected

            if accBleachIndex <= 0:
                # No acceptor signal detected either
                states.set_transition(n, 0, 
                        'donor_bleached', 0)
            
            elif accBleachIndex > 0:
                # The acceptor signal is terminated because of a donor bleach
                states.set_transition(n, accBleachIndex, 
                        'donor_bleached', bleachWindow)

        elif (donBleachIndex > 0) & (donBleachIndex < nrFrames):
            # A donor signal and bleach event were detected

            if donBleachIndex < accBleachIndex:
                # Something strange is going on - transition to a 
                # very high-fret state?
                states.set_transition(n, accBleachIndex, 
                        'donor_bleached', bleachWindow)

            elif donBleachIndex == accBleachIndex:
                # The acceptor signal is terminated because of a donor bleach
                states.set_transition(n, donBleachIndex, 
                        'donor_bleached', bleachWindow)

            elif donBleachIndex > accBleachIndex:
                # The usual situation
                states.set_transition(n, accBleachIndex, 
                        'acceptor_bleached', bleachWindow)
                states.set_transition(n, donBleachIndex, 
                        'donor_bleached', bleachWindow)

        elif donBleachIndex >= nrFrames:
            # The donor signal stayed throughout the measurement

            if accBleachIndex < nrFrames:
                # An acceptor bleach event was detected
                states.set_transition(n, accBleachIndex, 
                        'acceptor_bleached', bleachWindow)

            elif accBleachIndex >= nrFrames:
                # The acceptor signal also stayed throughout the measurement
                pass

        """
        #if not states.get(n).has('acceptor_present'):
        #if not accPresent:`
        #if accBleachIndex < donIllum[0].time.size:
        if accBleachIndex <= donBleachIndex:
            states.set_transition(n, accBleachIndex, 'acceptor_bleached', 3)
        else:
            states.set_transition(n, 0, 'acceptor_bleached', 3)
        
        #if (donBleachIndex > 0) & (donBleachIndex < donIllum[0].time.size):
        if donBleachIndex > 0:
            states.set_transition(n, donBleachIndex, 'donor_bleached', 3)
        """

        # == Direct acceptor illumination === #
        accData = accIllum[1].get(n).data
        nrAccDataPoints = accData.size

        if nrAccDataPoints > 0:
            accPresent = np.mean(accData) > traceSnr*np.std(accData)
        else:
            accPresent = False

        if accPresent:
            states.set_crit(n, acceptorLaserOn, 'acceptor_present')
        else:
            states.set_crit(n, acceptorLaserOn, 'misc')

        
    return states



def iteratively_get_bleach_index(data, detectionThreshold, zeroThreshold,
        bleachWindow, signalWindow):
    ## This first determines the bleach time using detectionThreshold.
    ## Then the data after the bleach is considered the baseline, but only
    ## if that baseline is close to zero as determined by zeroThresh
    data = np.array(data)

    if np.prod(data.shape)==data.size:
        data = np.reshape(data, data.size)

    originalData = np.array(data)

    originalIndex = get_bleach_index(data, detectionThreshold)
    currentIndex = originalIndex
    
    halfBleachWindow = np.round(bleachWindow/2.0)
    #if currentIndex < (data.size-BLEACH_WINDOW):
    #if currentIndex < data.size:
    for n in range(10):

        #print currentIndex
        if currentIndex < (data.size-(halfBleachWindow+signalWindow)):
            #baseline = np.mean(data[currentIndex+1:])
            intWin = int(currentIndex+halfBleachWindow)
            baseline = emean(data[intWin:intWin+signalWindow])
        else:
            baseline = emean(data)

        if baseline < zeroThreshold:
        #if baseline < detectionThreshold:
            data = data - baseline
            newIndex = get_bleach_index(data, detectionThreshold)
        else:
            break

        if newIndex > currentIndex:
            break

        if newIndex <= 0:
            startSignal = emean(originalData[:originalIndex])
            endSignal = emean(originalData[originalIndex:])

            if startSignal > zeroThreshold:
                if endSignal > zeroThreshold:
                    currentIndex = 2*data.size
                else:
                    currentIndex = originalIndex
            else:
                currentIndex = newIndex

            """
            if emean(originalData) < zeroThreshold:
                #currentIndex = newIndex
                currentIndex = originalIndex
            else:
                currentIndex = 2*data.size
            """
            break
         
        if np.abs(newIndex - currentIndex) < halfBleachWindow:
            #currentIndex = originalIndex
            break

        currentIndex = newIndex

    return currentIndex


def emean(data, *args):
    if data.size > 0:
        return np.mean(data, *args)
    else:
        return np.array([])

"""
def aget_signal_states(allTraces, descr,
        filterWindow,
        savGolOrder,
        traceSnr):
   
    #print('%s - Analyzing trace' % traceFile.filePath)

    laser532 = allTraces[descr.index('laser532')]
    laser638 = allTraces[descr.index('laser638')]
    
    traces = [allTraces[descr.index('donor')], 
            allTraces[descr.index('acceptor')]]

    donorLaserOn = laser532.data > 0
    acceptorLaserOn = laser638.data > 0
    
    #donor = traces[tracesHeader.index('donor')]
    #acceptor = traces[tracesHeader.index('acceptor')]

    donIllum = [trc.slice(donorLaserOn) for trc in traces]
    accIllum = [trc.slice(acceptorLaserOn) for trc in traces]

    # Smoothing polynomial filter used for analysis
    smoothFilter = SavGol(filterWindow, savGolOrder)
    #smoothFilter = MovingAverage(filterWindow)

    donIllumFilt = [trc.copy() for trc in donIllum]
    #misclib.list_apply(tracesFilt, smoothFilter.run)
    [smoothFilter.run(trc) for trc in donIllumFilt]
    
    # Determine noise RMS and signal threshold
    win = filterWindow

    #noiseRms = misclib.list_apply(donTraces, analysis.ts_get_rms_noise, win)
    noiseRms = [ts_get_rms_noise(trc, win) for trc in donIllum] 
    noiseThreshold = [nrms*traceSnr for nrms in noiseRms]
    
    signalWindow = get_signal_window(\
            *(donIllumFilt + noiseThreshold))
    accBleachWindow = get_acceptor_bleached_window(\
            *(donIllumFilt + noiseThreshold))
    #bleachWindow = get_all_bleached_window(\
    #        *(donIllumFilt + noiseThreshold))
    
    # Define states
    states = State(np.zeros(traces[0].data.shape))
    states.set_all_crit(np.ones(donIllum[0].data.shape), 'donor_bleached')
    states.set_all_crit(signalWindow, 'signal')
    states.set_all_crit(accBleachWindow, 'acceptor_bleached')
    
    # Introduce state transition zones
    nrTraces = states.get_nr_traces()
    bleachKernel = np.ones(3)
    for n in range(nrTraces):
        stt = states.get(n)
        accBleach = get_state_transitions(stt, ['signal', 'acceptor_bleached'])
        accBleach = np.convolve(accBleach, bleachKernel, mode='same')
        print accBleach.shape
        exit()



    # Acceptor presence
    accData = accIllum[1].data
    nrAccDataPoints = accData.shape[1]

    #states = State(traces[0])
    if nrAccDataPoints > 0:
        accPresent = np.mean(accData, axis=1) > traceSnr*np.std(accData, axis=1)
        fullAccPresent = accPresent
        for k in range(nrAccDataPoints-1):
            fullAccPresent = np.vstack((fullAccPresent, accPresent))
        fullAccPresent = np.transpose(fullAccPresent)
        
        crit = np.hstack((np.zeros(donIllum[0].data.shape), fullAccPresent))
        
        #states.data = np.hstack((stateCodeData, accPresent))
        
        states.set_all_crit(crit, 'acceptor_present')
    
    return states
"""


def get_bleach_index(data, thresh):
    if np.prod(data.shape)==data.size:
        data = np.reshape(data, data.size)

    sigHigh = data > thresh
    index = np.arange(data.size)
    indexHigh = index[sigHigh]
    farAway = 2*index[-1]
    if indexHigh.size == 0:
        return -farAway
    elif indexHigh[-1] == index[-1]:
        return index[-1] + farAway
    else:
        return indexHigh[-1]


"""
def ts_get_rms_noise(ts, win):
    nrTraces = ts.get_nr_traces()
    return np.array([np_get_rms_noise(ts.data[n,:], win) 
        for n in range(nrTraces)])

def np_get_rms_noise(data, win):
    return np.median(moving_rms(data, win))



def np_signal_is_high(data, thresh):
    return data > thresh

def np_get_bleach_time(data, time, thresh):
    sigH = np_signal_is_high(data, thresh)
    timeH = time[sigH]
    dt = np.mean(np.diff(time))
    if len(timeH) == 0:
        return 0
    elif timeH[-1] >= time[-1]:
        return time[-1]+dt
    else:
        return timeH[-1]

def np_has_bleached(data, time, thresh):
    bleachTime = np_get_bleach_time(data, time, thresh)
    if bleachTime > time[-1]:
        return time < 0
    else:
        return time >= bleachTime




def np_get_signal(don, acc, time, donThresh, accThresh):
    acceptorHasBleached = np_has_bleached(acc, time, accThresh)
    return acceptorHasBleached==False

def get_signal_window(don, acc, donThresh, accThresh):
    nrTraces = don.get_nr_traces()

    sig = []
    for n in range(nrTraces):
        signalWin = np_get_signal(don.get(n).data, acc.get(n).data, don.time, 
                donThresh[n], accThresh[n])
        sig.append(signalWin)
    #return TimeSeries(np.array(sig), don.time)
    return np.array(sig)




def np_get_acc_bleached(don, acc, time, donThresh, accThresh):
    acceptorHasBleached = np_has_bleached(acc, time, accThresh)
    donorHasBleached = np_has_bleached(don, time, donThresh)
    return np.logical_and(acceptorHasBleached, donorHasBleached==False)

def get_acceptor_bleached_window(don, acc, donThresh, accThresh):
    nrTraces = don.get_nr_traces()

    sig = []
    for n in range(nrTraces):
        signalWin = np_get_acc_bleached(don.get(n).data, 
                acc.get(n).data, don.time, 
                donThresh[n], accThresh[n])
        sig.append(signalWin)
    #return TimeSeries(np.array(sig), don.time)
    return np.array(sig)



def np_get_bleached(don, acc, time, donThresh, accThresh):
    donorHasBleached = np_has_bleached(don, time, donThresh)
    return donorHasBleached

def get_all_bleached_window(don, acc, donThresh, accThresh):
    nrTraces = don.get_nr_traces()

    sig = []
    for n in range(nrTraces):
        signalWin = np_get_bleached(don.get(n).data, 
                acc.get(n).data, don.time, 
                donThresh[n], accThresh[n])
        sig.append(signalWin)
    #return TimeSeries(np.array(sig), don.time)
    return np.array(sig)
"""
