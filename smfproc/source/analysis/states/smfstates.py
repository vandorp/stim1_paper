""" SMFSTATES """


# System packages
import os, shutil
import numpy as np

# Program modules
from smfproc.source.io.filetypes import TraceFile, SpecsFile
from smfproc.source.io.datatypes import State, SmfSelect
from smfproc.source.io.misc import get_output_filepaths as get_filepaths

# Program modules - local relative imports
from .analysis import get_trace_hmm, get_state_specs, get_state_sequence, get_gamma, \
        get_snr, get_state_transitions
from .adjustment import merge_states



### NB: Functions in this file can still use the TimeSeries class

stateDiff = lambda mu1, mu2: np.abs(mu1-mu2)/(mu1+mu2)

lxor = np.logical_xor
lor = np.logical_or
land = np.logical_and

def lloop(A, func):
    B = A[0]
    for n in range(1,len(A)):
        B = func(B, A[n])
    return B


def refine(args):
    """ Implements criteria for trace selection. """

    movieFilePath, coreSettings, traceSettings, stateSettings, selectSettings = args
    
    sourceFilePath, selectFilePath, specsFilePath, outputFilePath, snapFilePath = \
            get_filepaths(movieFilePath, coreSettings, 3)
    
    sourceFile = TraceFile(sourceFilePath)
    selectFile = SpecsFile(selectFilePath)
    specsFile = SpecsFile(specsFilePath)
    outputFile = TraceFile(outputFilePath)
    traces, header = sourceFile.read()
    selections = SmfSelect(selectFile.read())
    specs = specsFile.read()

    #selected = selections.get('selected')

    #autoRejected = selections.get(['no_fret', 'no_initial_signal'], np.logical_or)
    autoRejected = selections.get('no_fret')
    manualSelected = selections.get('manual')
    zeroFret = selections.get('zero_fret')
    toBeAnalyzed = np.logical_or(autoRejected==False, manualSelected)
    selected = autoRejected==False

    don, acc, fret, states = [traces[header.index(descr)] 
            for descr in ('donor', 'acceptor', 'fret', 'states')]
    states = State(states)


    nrTraces = don.get_nr_traces()
    totalTraces = np.sum(selected)
    cTrace = 0
    
    gRange = selectSettings['GAMMA_RANGE']
    gammas = np.zeros(nrTraces)
    snr = np.zeros(nrTraces)

    acceptorBleached = states.has('acceptor_bleach_event')
    donorBleached = states.has('donor_bleach_event')
    stableDonor = np.ones(nrTraces)==1
    validStateTransitions = np.ones(nrTraces)==1
    possibleSecondAcceptor = np.zeros(nrTraces)==1
    #possibleSecondDonor = np.zeros(nrTraces)==1
    
    #print('Trace\tSNR\tgamma\n')
    for n in range(nrTraces):
    #for n in [48]:
        if toBeAnalyzed[n]:
            cTrace += 1

            nStt = states.get(n)
            nDon = don.get(n)
            nAcc = acc.get(n)
            
            
            # 1 merge states with the bleach criterium
            # 2 determine the overall lowest state
            # 3 if it is below a criterium, and the trace ends 
            # with that state, then it is the bleached state


            # === Correct traces & states === #

                       
            
            
            # --- Relabel donor blink states ---

            donorStateCodes = nStt.get_codes('acceptor_bleached')
            #if acceptorBleached[n]:
            if len(donorStateCodes) > 1:

                donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                        donorStateCodes)
                donStd = get_state_specs(nDon.data, nStt.data, np.std, 
                        donorStateCodes, weighted=True)

                donorStates = \
                        nStt.data[nStt.get_all_crit('acceptor_bleached')]
                
                lowestStateMu = np.min(donMu)
                lowestStateCode = donorStateCodes[np.argmin(donMu)]

                #std = get_weighted_std(nDon, nStt, 'acceptor_bleached')
                if lowestStateMu < np.sum(donStd)*traceSettings['STATE_SNR']:
                    nStt.data[nStt.data==lowestStateCode] = nStt.codes['donor_blinked']

            
             # --- Correct donor bleach time ---
            # Get rid of 'blips'
            donorStateCodes = nStt.get_codes('acceptor_bleached')
            if len(donorStateCodes) > 1:
                if donorBleached[n]:
                    donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                            donorStateCodes)
                    donorStates = \
                            nStt.data[nStt.get_all_crit('acceptor_bleached')]
                    
                    blinkState = nStt.codes['donor_blinked']
                    bleachesFromBlinkState = donorStates[-1]==blinkState

                    if bleachesFromBlinkState:
                        donBleachEvents = get_state_transitions(nStt, 'acceptor_bleached')
                        newDonBleachIndex = np.where(donBleachEvents)[0][-2]
                        bleachWindow = traceSettings['BLEACH_WINDOW']

                        nStt.set_transition(0, newDonBleachIndex, 
                                'donor_bleached', bleachWindow)
                        states.set_transition(n, newDonBleachIndex, 
                                'donor_bleached', bleachWindow)

                    
            
            
            # --- Correct acceptor bleach time ---
            #correct_bleach_time() 
            
            # Determine signal-to-noise ratio
            sttCodes = nStt.get_codes('signal')
            sttSnr = get_snr(nDon.data+nAcc.data, nStt.data, sttCodes)
            
            #if (acceptorBleached[n] or zeroFret[n]):
            if zeroFret[n]:
                sttCodes = nStt.get_codes('acceptor_bleached')
                bleachSnr = get_snr(nDon.data+nAcc.data, nStt.data, sttCodes)
                #if bleachSnr < sttSnr:
                #    sttSnr = bleachSnr

                #if zeroFret[n]:
                sttSnr = bleachSnr

            snr[n] = sttSnr

            
            # --- Merge similar states ---
            try:
                states.data[n,:] = merge_states(nDon, nAcc, nStt, 
                        'signal', stateSettings['TRANSITION_SNR'], accountForSign=True)
            except:
                pass
            
            try:
                states.data[n,:] = merge_states(nDon, nAcc, nStt, 
                        'acceptor_bleached', stateSettings['TRANSITION_SNR'])
            except:
                pass
            
            # Correct baseline


            # === Refine selection === #
            # If there is only a single transition from an acceptor to
            # a lower fret state, then delete it (possible bleach of a 
            # second acceptor). Actually also do this for the donor, instead
            # of the below
            crit = nStt.get_all_crit('signal')
            sigData = nStt.data[crit]
            #if len(sigData) == 0:
            #    continue

            sttCodes, transTimes = get_state_sequence(sigData)
            accMu = get_state_specs(nAcc.data, nStt.data, np.mean, sttCodes)
            
            if len(sttCodes) > 1:
                singleTransitions = len(sttCodes) == len(list(set(sttCodes)))
                if singleTransitions:
                    decreasingFretStates = np.sum(np.diff(accMu)>=0)==0
                    if decreasingFretStates:
                        possibleSecondAcceptor[n] = True

            
            """
            crit = nStt.get_all_crit('acceptor_bleached')
            sttCodes, transTimes = get_state_sequence(nStt.data[crit])
            donMu = get_state_specs(nDon.data, nStt.data, np.mean, sttCodes)
            
            if len(sttCodes) == 2:
                
                singleTransitions = len(sttCodes) == len(list(set(sttCodes)))
                if singleTransitions:
                    decreasingFretStates = np.sum(diff(accMu)>=0)==0
                    if decreasingFretStates:
                        possibleSecondDonor = True
            
            elif len(sttCodes) > 2:
                possibleSecondDonor = True

            """

            # Only allow multiple donor states if the donor starts in the 
            # highest state and bleaches from that same highest state
            donorStateCodes = nStt.get_codes('acceptor_bleached')
            if len(donorStateCodes) == 2:
                if donorBleached[n]:
                    donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                            donorStateCodes)
                    donorStates = \
                            nStt.data[nStt.get_all_crit('acceptor_bleached')]

                    highestState = donorStateCodes[np.argmax(donMu)]
                    
                    startsInHighestState = donorStates[0]==highestState
                    bleachesFromHighestState = donorStates[-1]==highestState
                    
                    
                    if (startsInHighestState==False) or \
                            (bleachesFromHighestState==False):
                        stableDonor[n] = False



                    # ********************************************** #
                    # *** ALLOWING MILD DONOR FLUCTUATIONS HERE  *** #
                    # ********************************************** #
                    mildFluct = np.min(donMu) > 0.7*np.max(donMu)

                    if mildFluct:
                        stableDonor[n] = True
                        #pass



                else:
                    donorBlinked = nStt.has('donor_blinked')
                    if not donorBlinked:
                        stableDonor[n] = False

            elif len(donorStateCodes) > 2:
                stableDonor[n] = False

            
            """
            # Only allow multiple donor states if the donor starts in the 
            # highest state and there exists a transition from that highest 
            # state to a zero state
            donorStateCodes = nStt.get_codes('acceptor_bleached')
            if len(donorStateCodes) > 1:
                donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                        donorStateCodes)
                donorStates = \
                        nStt.data[nStt.get_all_crit('acceptor_bleached')]
                
                highestState = donorStateCodes[np.argmax(donMu)]
                blinkState = 
                
                startsInHighestState = donorStates[0]==highestState
                bleachedFromHighestState = \
                        donorBleached[n] and (donorStates[-1]==highestState)
                highestStateBlinks = True

                if len(donorStateCodes) > 3:
                    stableDonor[n] = False
                else:
                    if (startsInHighestState==False) or \
                            (bleachesFromHighestState==False):
                        stableDonor[n] = False

                    if highestStateBlinks:
                        stableDonor[n] = True

            """
            # Quit if there is no signal
            signalStates = \
                    nStt.data[nStt.get_all_crit('signal')]
            if signalStates.size < 3:
                validStateTransitions[n] = False
                continue

            stateCodes = nStt.get_codes('signal')
            nStates = len(stateCodes)
            gammaMatrix = np.zeros((nStates, nStates))
            for m in range(nStates):
                for k in range(nStates):
                    if not (m==k):
                        gammaMatrix[m,k] = get_gamma(nDon.data, nAcc.data, nStt.data, 
                                [stateCodes[m], stateCodes[k]])
            
            # Remove single-frame blips in the signal
            signalStates = \
                    nStt.data[nStt.get_all_crit('signal')]
            blipCrit = np.abs(np.diff(np.sign(np.diff(signalStates)))) > 1

            preBlipCrit = np.append(blipCrit, np.zeros(2)) == True
            blipCrit = np.append(0, np.append(blipCrit, 0)) == True
            
            preBlipStates = signalStates[preBlipCrit]
            blipStates = signalStates[blipCrit]
            for m in range(len(blipStates)):
                p = stateCodes.index(preBlipStates[m])
                k = stateCodes.index(blipStates[m])
                gamma = gammaMatrix[p,k]
                if gamma < 0:
                    crit = nStt.data[0,:] == blipStates[m]                    
                    nStt.data[0, crit] = preBlipStates[m]
                    states.data[n, crit] = preBlipStates[m]


            
            # Only allow state transitions within a specified range of gamma 
            # values
            stateCodes = nStt.get_codes('signal')
            nStates = len(stateCodes)
            gammaMatrix = np.zeros((nStates, nStates))
            for m in range(nStates):
                for k in range(nStates):
                    if m==k:
                        gammaMatrix[m,k] = np.nan
                    else:
                        gammaMatrix[m,k] = get_gamma(nDon.data, nAcc.data, nStt.data, 
                                [stateCodes[m], stateCodes[k]])

            """
            stateCodes = nStt.get_codes('signal')
            nStates = len(stateCodes)
            if nStates > 1:
                for m in range(nStates):
                    for k in range(nStates):
                        if not (m==k):
                            gamma = get_gamma(nDon.data, nAcc.data, nStt.data, 
                                    [stateCodes[m], stateCodes[k]])
                            
                            #if (gamma < gRange[0]) or (gamma > gRange[1]):
                            if (gamma < 0) or (gamma > 2*gRange[1]):
                                validStateTransitions[n] = False
            """
            M = np.array(gammaMatrix)
            M[np.isnan(M)] = np.mean(gRange)
            if (np.sum(M < 0.5*gRange[0]) > 0) or \
                    (np.sum(M > 2*gRange[1]) > 0):
                validStateTransitions[n] = False


            
            # Determine bleach gamma values
            if acceptorBleached[n]:
                signalStates = \
                            nStt.data[nStt.get_all_crit('signal')]
                bleachedStates = \
                            nStt.data[nStt.get_all_crit('acceptor_bleached')]

                signalStateCodes = nStt.get_codes('signal')
                accMu = get_state_specs(nAcc.data, nStt.data, np.mean, 
                        signalStateCodes)
                maxFretState = signalStateCodes[np.argmax(accMu)]


                donorStateCodes = nStt.get_codes('acceptor_bleached')
                donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                        donorStateCodes)
                
                if (len(bleachedStates) > 0): # and (len(signalStates) > 0):


                    # Calculate gamma based on signals just before and after
                    # the acceptor bleach event.
                    #gamma = get_gamma(nDon.data, nAcc.data, nStt.data,
                    #        [signalStates[-1], bleachedStates[0]])

                    # Calculate gamma based on signal in the highest FRET
                    # state and signal just after the acceptor bleach event.
                    #gamma = get_gamma(nDon.data, nAcc.data, nStt.data,
                    #        [maxFretState, bleachedStates[0]])

                    # Calculate gamma based on signal in the highest FRET
                    # state and highest donor signal after the acceptor bleach event.

                    if len(donMu) > 0:
                        maxDonorState = donorStateCodes[np.argmax(donMu)]
                    else:
                        maxDonorState = bleachedStates[0]

                    gamma = get_gamma(nDon.data, nAcc.data, nStt.data,
                            [maxFretState, maxDonorState])

                    gammas[n] = gamma

                else:
                    # The bleached state is too short to measure

                    if gammaMatrix.size > 1:
                        # if there are transitions in the trace, then calculate gamma
                        gammas[n] = np.nanmean(gammaMatrix)
                    else:
                        # if no transitions, then use signal contained within the bleach window
                        #accBleachedState = nStt.data[nStt.get_all_crit('acceptor_bleach_event')]

                        accMu = get_state_specs(nAcc.data, nStt.data, np.mean, 
                                [maxFretState])
                        donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                                [maxFretState])

                        highFret = accMu[0] > donMu[0]

                        if highFret:
                            gamma = get_gamma(nDon.data, nAcc.data, nStt.data,
                                    [maxFretState, 2])
                            gammas[n] = gamma




            # If a trace has no acceptor bleach, but it does have state
            # transitions to near zero, include it and calculate gamma
            if not acceptorBleached[n]:
                stateCodes = nStt.get_codes('signal')
                if len(stateCodes) > 1:
                    donMu = get_state_specs(nDon.data, nStt.data, np.mean, 
                            stateCodes)
                    accMu = get_state_specs(nAcc.data, nStt.data, np.mean, 
                            stateCodes)

                    minState = np.argmin(accMu)
                    maxState = np.argmax(accMu)
                    
                    condition1 = accMu[minState] < 0.25*accMu[maxState]
                    condition2 = donorBleached[n]

                    if condition1 or condition2:
                        acceptorBleached[n] = True
                        gamma = get_gamma(nDon.data, nAcc.data, nStt.data,
                                [stateCodes[maxState], stateCodes[minState]])
                        gammas[n] = gamma
            

            # Set a threshold for acceptor presence
            if zeroFret[n]:
                stateCodes = nStt.get_codes('acceptor_present')
                accMu = get_state_specs(nAcc.data, nStt.data, np.mean, 
                            stateCodes)

                stateCodes = nStt.get_codes('acceptor_bleached')
                accStd = get_state_specs(nAcc.data, nStt.data, np.std, 
                            stateCodes)
                accN = get_state_specs(nAcc.data, nStt.data, len, 
                            stateCodes)
                
                if (len(accMu) > 0) and (len(accStd) > 0):
                    accStd = np.sum(accStd*accN)/np.sum(accN)

                    if accMu < 3*accStd:
                        zeroFret[n] = False
                else:
                    zeroFret[n] = False


            signalStateCodes = nStt.get_codes('signal')
            #if acceptorBleached[n]:
            if len(signalStateCodes) > 1:

                accMu = get_state_specs(nAcc.data, nStt.data, np.mean, 
                        signalStateCodes)
                accStd = get_state_specs(nAcc.data, nStt.data, np.std, 
                        signalStateCodes, weighted=True)

                accStates = \
                        nStt.data[nStt.get_all_crit('signal')]
                
                lowestStateMu = np.min(accMu)
                lowestStateCode = signalStateCodes[np.argmin(accMu)]

                #std = get_weighted_std(nDon, nStt, 'acceptor_bleached')
                if lowestStateMu < np.sum(accStd)*traceSettings['STATE_SNR']:

                    zeroData = nStt.data[nStt.data==lowestStateCode]
                    if len(zeroData) > 0.5*len(accStates):
                        acceptorBleached[n] = False
                        zeroFret[n] = True
            

        cIndex = n
        cSnr = snr[n]
        cGamma = gammas[n]
        #print('%d\t%f\t%f'%(cIndex, cSnr, cGamma))



    
    # === write new states to file === #
    outputFile.write_block(states, 'states')
    shutil.copyfile(outputFilePath, snapFilePath)
    
    
    # === Store parameters === #
    specs.set(gammas, 'gamma')
    specs.set(snr, 'snr')
    specsFile.write(specs)
    
    
    # === Make selections === #
    

    
    # >>> === SELECTION CRITERIA === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #newRejections = SmfSelect()

    #rejectDescr = ['no_fret']
    
    """
    if selectSettings['REQUIRE_ACCEPTOR_BLEACH']:
        crit = land(selected, acceptorBleached==False)
        newRejections.set(crit, 'no_acceptor_bleach')
    
    if selectSettings['REQUIRE_DONOR_BLEACH']:
        crit = land(selected, donorBleached==False)
        newRejections.set(crit, 'no_donor_bleach')
    """
    
    # If the gamma could not be determined, the trace is _not_ accepted
    G = np.array(gammas)
    G[np.isnan(G)] = gRange[0]-1
    validAcceptorBleach = land((G >= gRange[0]), (G <= gRange[1]))
    

    # if the snr could not be determined, the trace is _not_ accepted
    snr[np.isnan(snr)] = selectSettings['MINIMUM_SNR']-1
    validSnr = np.logical_and(snr >= selectSettings['MINIMUM_SNR'], snr < 15)
    #validSnr = np.logical_or(tempSnr >= selectSettings['MINIMUM_SNR'], np.isnan(snr))

    okaySnr = np.logical_and(snr >= 15, donorBleached)

    validSnr = np.logical_or(validSnr, okaySnr)

    """
    newRejections.set(land(selected,stableDonor==False), 'rej_unstable_donor')
    #newRejections.set(land(selected,possibleSecondDonor), 'second_donor')
    newRejections.set(land(selected,possibleSecondAcceptor), 'rej_second_acceptor')
    newRejections.set(land(selected,validStateTransitions==False), 'rej_invalid_transitions')
    newRejections.set(land(selected,validAcceptorBleach==False), 'rej_invalid_acceptor_bleach')
    newRejections.set(land(selected,validSnr==False), 'rej_bad_snr')

    rejectDescr.extend(newRejections.descr)
    selections.merge(newRejections)
    

    rejected = selections.get(rejectDescr, lor)
    """

    noFret = selections.get('no_fret')

    selections.set(land(selected,stableDonor==False), 'rej_unstable_donor')
    #newRejections.set(land(selected,possibleSecondDonor), 'second_donor')
    selections.set(land(selected,possibleSecondAcceptor), 'rej_second_acceptor')
    selections.set(land(selected,validStateTransitions==False), 'rej_invalid_transitions')
    selections.set(land(selected,validAcceptorBleach==False), 'rej_invalid_acceptor_bleach')
    selections.set(land(selected,validSnr==False), 'rej_bad_snr')


    # >>> === CATEGORIES === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    crit = lloop((acceptorBleached==False, donorBleached==False, noFret==False,
        stableDonor==True, possibleSecondAcceptor==False, zeroFret==False,
        validStateTransitions==True, validSnr==True), land)
    selections.set(crit, 'cat_no_bleach')

    crit = lloop((acceptorBleached==False, donorBleached==True, noFret==False,
        stableDonor==True, possibleSecondAcceptor==False, zeroFret==False,
        validStateTransitions==True, validSnr==True), land)
    selections.set(crit, 'cat_donor_bleach_only')

    crit = lloop((acceptorBleached==True, donorBleached==False, noFret==False,
        stableDonor==True, possibleSecondAcceptor==False, 
        validStateTransitions==True, validSnr==True, 
        validAcceptorBleach==True), land)
    selections.set(crit, 'cat_acceptor_bleach_only')

    crit = lloop((acceptorBleached==True, donorBleached==True, noFret==False,
        stableDonor==True, possibleSecondAcceptor==False, 
        validStateTransitions==True, validSnr==True, 
        validAcceptorBleach==True), land)
    selections.set(crit, 'cat_both_bleached')
    
    crit = lloop((donorBleached==False, zeroFret==True, noFret==False,
        stableDonor==True, validSnr==True), land)
    selections.set(crit, 'cat_zero_fret_no_bleach') 

    crit = lloop((donorBleached==True, zeroFret==True, noFret==False,
        stableDonor==True, validSnr==True), land)
    selections.set(crit, 'cat_zero_fret_donor_bleach') 
    
    # Make final selections
    #zeroFretSelected = land(land(land(zeroFret, stableDonor), validSnr), donorBleached)
    #selections.set(zeroFretSelected, 'selected_zero_fret')


    #rejected = selections.get(rejectDescr, lor)
    #crit = land(lor(rejected==False, zeroFretSelected), manualSelected==False)
    crit = selections.get([
        'cat_no_bleach',
        'cat_donor_bleach_only',
        'cat_acceptor_bleach_only', 
        'cat_both_bleached',
        'cat_zero_fret_no_bleach',
        'cat_zero_fret_donor_bleach'], lor)
    selections.set(crit, 'auto')

    finalSelected = selections.get(['manual', 'auto'], lxor)
    selections.set(finalSelected, 'selected')

    # Write to file
    selectFile.write(selections)

#from tqdm import tqdm
#import progressbar


def analyze(args):
    """ Top level analysis.
    
    Handles I/O, some data prep.
    """

    # --- GET DATA ---
    movieFilePath, coreSettings, traceSettings = args
    sourceFilePath, selectFilePath, specsFilePath, outputFilePath, snapFilePath = \
            get_filepaths(movieFilePath, coreSettings, 2)

    sourceFile = TraceFile(sourceFilePath)
    selectFile = SpecsFile(selectFilePath)
    outputFile = TraceFile(outputFilePath)

    traces, header = sourceFile.read()
    
    selections = SmfSelect(selectFile.read())
    #autoRejected = selections.get(['no_fret', 'no_initial_signal'], np.logical_or)
    autoRejected = selections.get('no_fret')
    selected = autoRejected==False

    don, acc, fret, states = [traces[header.index(descr)] 
            for descr in ('donor', 'acceptor', 'fret', 'states')]
    states = State(states)
    # -------------

    #movieName = os.path.basename(movieFilePath)
    movieName = os.sep.join(movieFilePath.split(os.sep)[-3:])
    print('Fitting HMM states for movie %s... '%movieName)
    
    nrTraces = don.get_nr_traces()
    totalTraces = np.sum(selected)
    cTrace = 0
    
    #progress = progressbar.ProgressBar()
    #progress = progressbar.ProgressBar(widgets=[
    #    progressbar.Bar('=', '[', ']'), ' ',
    #    progressbar.Percentage(), ' ',
    #    progressbar.ETA()])

    #for n in progress(range(nrTraces)):
    for n in range(nrTraces):
        if selected[n]:
            cTrace += 1
            #print 'Processing trace #%d (%d/%d)'%(n, cTrace, totalTraces)

            
            newStates, params = get_trace_hmm([acc.get(n), don.get(n)], 
                    states.get(n), 'signal', traceSettings['MAX_NR_STATES'])
            states.data[n,:] = newStates.data[0,:]
            
            newStates, params = get_trace_hmm([acc.get(n), don.get(n)], 
                    states.get(n), 'acceptor_bleached', traceSettings['MAX_NR_STATES'])
            states.data[n,:] = newStates.data[0,:]

    # write new states to file
    outputFile.write_block(states, 'states')
    shutil.copyfile(outputFilePath, snapFilePath)

