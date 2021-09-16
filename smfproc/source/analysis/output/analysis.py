""" Analysis output and figures. """

# System packages
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# Program modules
from smfproc.source.io.filetypes import SmfSelect
from smfproc.source.io.datatypes import State
from smfproc.source.analysis.states.analysis import get_state_data, get_state_specs, get_state_model
from smfproc.source.analysis.traces.analysis import correct_fret




""" SETTINGS """

HISTOGRAM_RANGE = (-0.25, 1.25)
HISTOGRAM_NR_BINS = 40
HISTOGRAM_COLOR = 'white'


def draw_fret_histogram(ax, bins, cnt, *args, **kwargs):
    width = np.mean(np.diff(bins))
    h = ax.bar(bins, cnt, width=width, color=HISTOGRAM_COLOR, align='center', 
            *args, **kwargs)

    ax = plt.gca()
    ax.set_xlim(HISTOGRAM_RANGE)
    ax.set_xlabel('FRET ratio')
    ax.set_ylabel('Count (norm)')
    return h


def get_fret_histograms(traces, norm_indiv=False, norm_total=False, nBins=HISTOGRAM_NR_BINS):

    allCnt = np.zeros(nBins)
    bins = np.linspace(HISTOGRAM_RANGE[0], HISTOGRAM_RANGE[1], nBins)

    if len(traces) == 0:
        return allCnt, bins, 1.0
    
    #allCnt = []

    for trc in traces:

        if len(trc) == 0:
            continue

        combinedStates = np.hstack([stt.data for stt in trc])

        cnt, binEdges = np.histogram(combinedStates,
                range=HISTOGRAM_RANGE, 
                bins=nBins)
        
        if norm_indiv:
            if np.sum(cnt) > 0:
                cnt = (1.0*cnt)/np.sum(cnt)

        #allCnt.append(cnt)
        allCnt = allCnt + cnt
    
    #allCnt = np.sum(np.array(allCnt), axis=0)
    
    normCoeff = np.sum(allCnt)
    if (norm_total and (normCoeff > 0)):
        allCnt = (1.0*allCnt)/normCoeff
    
    #bins = (binEdges[:-1]+binEdges[1:])/2
    return allCnt, bins, normCoeff


"""
def save_analysis_data(data, outDir, flag, header=''):
    dataOutPath = os.path.join(outDir, flag+'.hist')
    #header = 'bins fret fret_norm fret_corr fret_norm_corr'
    np.savetxt(dataOutPath, data, '%.3f', header=header)
"""

def save_analysis_figure(fig, outPath):
    #figOutPath = os.path.join(outDir, flag+'.png')
    plt.figure(fig)
    plt.savefig(outPath)
    plt.close(fig)


def listsubset(A, c):
    return [A[n] for n in range(len(A)) if c[n]]


def get_categories_histograms_data(fret, states, selections, stype='auto'):
    signal = get_state_data(fret, states, descr='signal')

    selDescr = [
            'cat_no_bleach',
            'cat_donor_bleach_only',
            'cat_acceptor_bleach_only',
            'cat_both_bleached']
    
    cnt = []
    descr = []
    for sd in selDescr:
        crit = selections.get(sd)
        if stype == 'manual':
            crit_manual = selections.get('manual')
            crit = np.logical_and(crit, crit_manual)
        trc = fret.select(crit)
        stt = states.select(crit)
        signal = get_state_data(trc, stt, descr='signal')

        #cntFret, bins, _ = get_fret_histograms(
        #        listsubset(signal, selections.get(sd)))
        cntFret, bins, _ = get_fret_histograms(signal)
        cntFretNorm, bins, _ = get_fret_histograms(signal, norm_indiv=True)

        cnt.extend([cntFret, cntFretNorm])
        descr.extend([sd, sd+'_norm'])

    """
    for sd in selDescr:
        cntFretNorm, bins, _ = get_fret_histograms(
                listsubset(signal, selections.get(sd)), norm_indiv=True)
        
        
        cnt.append(cntFretNorm)
        descr.append(sd+'_norm')
    """

    return np.array(cnt), bins, descr


def maxcomp(data, ref):
    dmax = np.max(data)
    if dmax > ref:
        return dmax
    else:
        return ref


def make_allpoints_histogram_figure(ax, cntFret, bins, descr):
    nrPoints = np.sum(cntFret)
    
    if nrPoints > 0:
        cntFret = 1.0*cntFret/nrPoints

    h = draw_fret_histogram(ax, bins, cntFret)
    #ymax = maxcomp(cntFret/nrPoints, ymax)
    
    ax.set_ylabel('Count (norm)')
    ax.set_xlabel('')
    ax.set_xlim(HISTOGRAM_RANGE)
    ax.set_xticklabels([])
    ax.set_title('%s (N=%d)'%(descr, round(nrPoints)), fontsize=10)
    return h


def make_norm_histogram_figure(ax, cntFret, bins, descr):
    h = make_allpoints_histogram_figure(ax, cntFret, bins, descr)
    ax.set_ylabel('')
    ax.set_yticklabels([])
    return h


def get_categories_histograms_figure(fig, cnt, bins, descr):
    #fig = plt.figure(1, figsize=(6,9))
    ymax = 0

    cnt = np.transpose(cnt)
    
    # === NO DONOR AND NO ACCEPTOR BLEACH === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cntFret = cnt[:,0]
    cntFretNorm = cnt[:,1]
    
    ax = fig.add_subplot(4,2,1)
    h = make_allpoints_histogram_figure(ax, cntFret, bins, descr[0][4:])

    ax = fig.add_subplot(4,2,2)
    h = make_norm_histogram_figure(ax, cntFretNorm, bins, descr[1][4:])
    
    # === DONOR BLEACH ONLY === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cntFret = cnt[:,2]
    cntFretNorm = cnt[:,3]

    ax = fig.add_subplot(4,2,3)
    h = make_allpoints_histogram_figure(ax, cntFret, bins, descr[2][4:])

    ax = fig.add_subplot(4,2,4)
    h = make_norm_histogram_figure(ax, cntFretNorm, bins, descr[3][4:])

    # === ACCEPTOR BLEACH ONLY === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cntFret = cnt[:,4]
    cntFretNorm = cnt[:,5]
    
    ax = fig.add_subplot(4,2,5)
    h = make_allpoints_histogram_figure(ax, cntFret, bins, descr[4][4:])

    ax = fig.add_subplot(4,2,6)
    h = make_norm_histogram_figure(ax, cntFretNorm, bins, descr[5][4:])

    # === DONOR AND ACCEPTOR BLEACH === <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cntFret = cnt[:,6]
    cntFretNorm = cnt[:,7]

    ax = fig.add_subplot(4,2,7)
    h = make_allpoints_histogram_figure(ax, cntFret, bins, descr[6][4:])
    ax.set_xlabel('FRET ratio')
    ax.set_xticklabels(ax.get_xticks(), fontsize=10)

    ax = fig.add_subplot(4,2,8)
    h = make_norm_histogram_figure(ax, cntFretNorm, bins, descr[7][4:])
    ax.set_xlabel('FRET ratio')
    ax.set_xticklabels(ax.get_xticks(), fontsize=10)

    ymax = max([a.get_ylim()[1] for a in fig.axes])
    [a.set_ylim([0, 1.2*ymax]) for a in fig.axes]
    #return 1



def stack_fret_traces(donor, acceptor, states, selections, gamma):
    crit1 = selections.get('cat_acceptor_bleach_only')
    crit2 = selections.get('cat_both_bleached')
    crit = np.logical_or(crit1, crit2)

    # Fluctuating traces only
    nrTraces = states.get_nr_traces()
    nrStates = np.array([len(states.get(n).get_codes('signal')) 
            for n in range(nrTraces)])
    crit3 = nrStates > 1
    #crit = np.logical_and(crit, crit3)

    donorSelected = donor.select(crit)
    acceptorSelected = acceptor.select(crit)
    statesSelected = states.select(crit)
    gamma = gamma[crit]
    
    donorCorr = donorSelected.copy()
    for m in range(donorCorr.get_nr_traces()):
        donorCorr.data[m,:] = donorCorr.data[m,:]*gamma[m]
    
    crit = statesSelected.get_all_crit('signal')

    #don = get_state_data(donorCorr, statesSelected, descr='signal')
    #acc = get_state_data(acceptorSelected, statesSelected, descr='signal')
    
    duration = 30
    allData = []
    for m in range(donorCorr.get_nr_traces()):
        don = donorCorr.data[m, crit[m,:]]
        acc = acceptorSelected.data[m, crit[m,:]]
        index = (m+1)*np.ones(don.size)
        data = np.transpose(np.vstack((index, don, acc)))
        if don.size >= duration:
            allData.append(data)
    if len(allData)>0:
        data = np.vstack(allData)
    else:
        data = np.array([])

    return data


def get_hmm_histograms_data(fret, states, selections, gamma, stype):
    if stype == 'manual':
        crit = selections.get('manual') > 0
    else:

        crit1 = selections.get('cat_acceptor_bleach_only')
        crit2 = selections.get('cat_both_bleached')
        crit = np.logical_or(crit1, crit2)

    fretSelected = fret.select(crit)
    statesSelected = states.select(crit)
    gamma = gamma[crit]
    
    fretCorr = fretSelected.copy()
    correct_fret(fretCorr.data, gamma)
    
    model = []
    modelCorr = []

    nrTraces = fretCorr.get_nr_traces()
    for index in range(nrTraces):
        fretModel = get_state_model(
                fretSelected.get(index), statesSelected.get(index), descr='signal')
        fretCorrModel = get_state_model(
                fretCorr.get(index), statesSelected.get(index), descr='signal')
        
        model.append([fretModel])
        modelCorr.append([fretCorrModel])

    cntFret, bins, _ = get_fret_histograms(model)
    cntCorr, bins, _ = get_fret_histograms(modelCorr)
    
    # === Combine histograms
    cnt = np.vstack((cntFret, cntCorr))
    descr = ['model', 'model_corr']
    return cnt, bins, descr



def get_fret_histograms_data(fret, states, selections, gamma, stype='auto'):
    if stype == 'manual':
        crit = selections.get('manual') > 0
    else:

        crit1 = selections.get('cat_acceptor_bleach_only')
        crit2 = selections.get('cat_both_bleached')
        crit = np.logical_or(crit1, crit2)

    fretSelected = fret.select(crit)
    statesSelected = states.select(crit)
    gamma_ = gamma[crit]
    
    fretCorr = fretSelected.copy()
    correct_fret(fretCorr.data, gamma_)

    signal = get_state_data(fretSelected, statesSelected, descr='signal')
    signalCorr = get_state_data(fretCorr, statesSelected, descr='signal')

    cntFret, bins, _ = get_fret_histograms(signal)
    cntNorm, bins, indivNormCoeff = get_fret_histograms(signal, norm_indiv=True)
    cntCorr, bins, _ = get_fret_histograms(signalCorr)
    cntCorrNorm, bins, _ = get_fret_histograms(signalCorr, norm_indiv=True)
    
    """
    # === Combine all categories
    crit3 = selections.get('cat_donor_bleach_only')
    crit4 = selections.get('cat_no_bleach')
    crit = (crit + crit3 + crit4) > 0

    fretSelected = fret.select(crit)
    statesSelected = states.select(crit)
    
    signal = get_state_data(fretSelected, statesSelected, descr='signal')
    cntAllNorm, bins, indivNormCoeff = get_fret_histograms(signal, 
            norm_indiv=True, norm_total=False)
    """

    # === Get histogram of zero-fret traces
    crit = np.logical_or(selections.get('cat_zero_fret_donor_bleach'),
            selections.get('cat_zero_fret_no_bleach'))
    if stype == 'manual':
        crit_manual = selections.get('manual')
        crit = np.logical_and(crit, crit_manual)


    #crit = selections.get('cat_zero_fret_donor_bleach')
    fretZero = fret.select(crit)
    statesZero = states.select(crit)

    zeroFretSignal = get_state_data(fretZero, statesZero, descr='acceptor_bleached')
    cntZeroFretNorm, _, _ = get_fret_histograms(zeroFretSignal, 
            norm_indiv=True, norm_total=False)

    if len(cntZeroFretNorm) == 0:
        cntZeroFretNorm = np.zeros(HISTOGRAM_NR_BINS)
    
    """
    # === Get manual selected traces
    crit = selections.get('manual') > 0

    fretSelected = fret.select(crit)
    statesSelected = states.select(crit)
    gamma_ = gamma[crit]
    
    fretCorr = fretSelected.copy()
    correct_fret(fretCorr.data, gamma_)

    signal = get_state_data(fretSelected, statesSelected, descr='signal')
    signalCorr = get_state_data(fretCorr, statesSelected, descr='signal')

    cntFret_manual, bins, _ = get_fret_histograms(signal)
    cntNorm_manual, bins, indivNormCoeff = get_fret_histograms(signal, norm_indiv=True)
    cntCorr_manual, bins, _ = get_fret_histograms(signalCorr)
    cntCorrNorm_manual, bins, _ = get_fret_histograms(signalCorr, norm_indiv=True)
    """

    # === Combine histograms
    #cnt = np.vstack((cntFret, cntCorr, cntNorm, cntCorrNorm, cntAllNorm, cntZeroFretNorm, cntFret_manual, cntNorm_manual, cntCorr_manual, cntCorrNorm_manual))
    #descr = ['fret', 'fret_corr', 'fret_norm', 'fret_norm_corr', 'fret_all_norm', 'zero_fret', 'manual_fret', 'manual_fret_norm', 'manual_fret_corr', 'manual_fret_corr_norm']
    cnt = np.vstack((cntFret, cntCorr, cntNorm, cntCorrNorm, cntZeroFretNorm))
    descr = ['fret', 'fret_corr', 'fret_norm', 'fret_norm_corr', 'zero_fret']

    return cnt, bins, descr


def norm_hist(y):
    if np.sum(y)>0:
        return y/np.sum(y)
    else:
        return y


def get_hmm_histograms_figures(ax1, cnt, bins, descr, fontSize=12):
    cnt[np.isnan(cnt)] = 0
    
    model = cnt[:,0]
    modelCorr = cnt[:,1]
    
    # Sample down
    bins = bins[::2]
    modelCorr = modelCorr[::2] + modelCorr[1::2]


    draw_fret_histogram(ax1, bins, norm_hist(modelCorr))
    #ax1.plot(bins, norm_hist(modelCorr), '-k')
    #ax1.plot(bins, norm_hist(model), '-r')
    ax1.set_xlabel('FRET ratio')
    ax1.set_ylabel('Count')
    ax1.set_xlim(HISTOGRAM_RANGE)
    ax1.legend({'N = %d'%np.sum(model)}, fontsize=fontSize)


def get_fret_histograms_figures(ax1, ax2, cnt, bins, descr, fontSize=12):
    #cnt = np.transpose(cnt)

    cnt[np.isnan(cnt)] = 0
    
    fret = cnt[:,0]
    fretNorm = cnt[:,2]
    fretCorr = cnt[:,1]
    fretCorrNorm = cnt[:,3]
    zeroFret = cnt[:,4]
    
    zeroNorm = np.sum(cnt[:,1]+cnt[:,4])

    #fig2 = plt.figure(2)
    #ax = plt.gca()
    draw_fret_histogram(ax1, bins, norm_hist(fretCorr))
    ax1.plot(bins, norm_hist(fretCorr), '-k')
    ax1.plot(bins, norm_hist(fret), '-r')
    ax1.set_xlabel('FRET ratio')
    ax1.set_ylabel('Count')
    ax1.set_xlim(HISTOGRAM_RANGE)
    ax1.legend({'N = %d'%np.sum(fret)}, fontsize=fontSize)
    #ax1.set_title('N=%d'%np.sum(fret), fontsize=12)
    
    
    #fig3 = plt.figure(3)
    #ax = plt.gca()
    draw_fret_histogram(ax2, bins, norm_hist(fretCorrNorm))
    #plot_stacked_bars(ax2, bins, np.vstack((fretCorrNorm, zeroFret)), stacked=True, colors=('w', [0.8,0.8,0.8]))
    ax2.plot(bins, norm_hist(fretCorrNorm), '-k')
    ax2.plot(bins, norm_hist(fretNorm), '-r')
    ax2.set_xlabel('FRET ratio')
    ax2.set_ylabel('Count (norm)')
    ax2.set_xlim(HISTOGRAM_RANGE)
    ax2.legend({'N = %d'%np.sum(fretCorrNorm)}, fontsize=fontSize)
    #ax2.set_title('$N_{signal}=%d$, $N_{zero}=%d$'%(np.sum(fretCorrNorm), np.sum(zeroFret)), fontsize=10)
    
    """
    #draw_fret_histogram(ax3, bins, fretCorrNorm)
    #plot_stacked_bars(ax3, bins, np.vstack((fretAll/zeroNorm, zeroFret/zeroNorm)), stacked=True, colors=('w', [0.8,0.8,0.8]))
    plot_stacked_bars(ax3, bins, np.vstack((norm_hist(fretCorr), norm_hist(zeroFret))), stacked=True, colors=('w', [0.8,0.8,0.8]))

    ax3.set_xlabel('FRET ratio')
    ax3.set_ylabel('Count')
    ax3.set_xlim(HISTOGRAM_RANGE)
    #ax3.legend(['N=%d'%np.sum(fretAll), 'N=%d'%np.sum(zeroFret)], fontsize=fontSize)
    """

"""
def get_fret_histograms_figures_manual(ax1, ax2, cnt, bins, descr, fontSize=12):
    #cnt = np.transpose(cnt)

    cnt[np.isnan(cnt)] = 0

    fret = cnt[:,6]
    fretNorm = cnt[:,7]
    fretCorr = cnt[:,8]
    fretCorrNorm = cnt[:,9]

    #fig2 = plt.figure(2)
    #ax = plt.gca()
    draw_fret_histogram(ax1, bins, norm_hist(fretCorr))
    ax1.plot(bins, norm_hist(fretCorr), '-k')
    ax1.plot(bins, norm_hist(fret), '-r')
    ax1.set_xlabel('FRET ratio')
    ax1.set_ylabel('Count')
    ax1.set_xlim(HISTOGRAM_RANGE)
    ax1.legend({'N = %d'%np.sum(fret)}, fontsize=fontSize)
    #ax1.set_title('N=%d'%np.sum(fret), fontsize=12)
    
    
    #fig3 = plt.figure(3)
    #ax = plt.gca()
    draw_fret_histogram(ax2, bins, norm_hist(fretCorrNorm))
    #plot_stacked_bars(ax2, bins, np.vstack((fretCorrNorm, zeroFret)), stacked=True, colors=('w', [0.8,0.8,0.8]))
    ax2.plot(bins, norm_hist(fretCorrNorm), '-k')
    ax2.plot(bins, norm_hist(fretNorm), '-r')
    ax2.set_xlabel('FRET ratio')
    ax2.set_ylabel('Count (norm)')
    ax2.set_xlim(HISTOGRAM_RANGE)
    ax2.legend({'N = %d'%np.sum(fretCorrNorm)}, fontsize=fontSize)
    #ax2.set_title('$N_{signal}=%d$, $N_{zero}=%d$'%(np.sum(fretCorrNorm), np.sum(zeroFret)), fontsize=10)
"""

def get_nrstates_histograms_data(fret, states, selections, gamma, stype='auto'):
    if stype == 'manual':
        crit = selections.get('manual') > 0
    else:

        crit1 = selections.get('cat_acceptor_bleach_only')
        crit2 = selections.get('cat_both_bleached')
        crit = np.logical_or(crit1, crit2)

    fretSelected = fret.select(crit)
    statesSelected = states.select(crit)
    gamma = gamma[crit]
    
    fretCorr = fretSelected.copy()
    correct_fret(fretCorr.data, gamma)

    nrTraces = statesSelected.get_nr_traces()
    nrStates = np.array([len(statesSelected.get(n).get_codes('signal')) 
            for n in range(nrTraces)])
    
    dt = np.mean(np.diff(fret.time))
    duration = np.sum(statesSelected.get_all_crit('signal'), axis=1)*dt

    allStack = []
    for n in range(1,5):
        crit = np.logical_and(nrStates==n, duration>=30)
        signal = get_state_data(
                fretCorr.select(crit),
                statesSelected.select(crit),
                descr='signal')
        cntFret, bins, _ = get_fret_histograms(signal, norm_indiv=True)
        allStack.append(cntFret)
    
    allStack = np.array(allStack)
    
    descr = ['one_state', 'two_states', 'three_states', 'four_states']

    return allStack, bins, descr

   
def get_nrstates_histograms_figure(ax, cnt, bins, descr, fontSize=12):
    cnt = np.transpose(cnt)
    #cnt[np.isnan(cnt)] = 0

    nrsTraces = np.array(np.round(np.sum(cnt, axis=1)), dtype=int)
    cnt = np.flipud(cnt)
    nrTraces = np.round(np.sum(cnt))
    
    if nrTraces > 0:
        cnt = cnt/nrTraces

    #fig4 = plt.figure(4)
    #ax = plt.gca()
    #plot_stacked_bars(ax, bins, cnt, stacked=True, colors=('m','b','g','r'))
    plot_stacked_bars(ax, bins, cnt, stacked=True, colors=(
        [0.4,0.4,0.4],
        [0.6,0.6,0.6],
        [0.8,0.8,0.8],
        'w'))

    ax.set_xlabel('FRET ratio')
    ax.set_ylabel('Count (norm)')
    ax.set_xlim(HISTOGRAM_RANGE)
    #ax.set_title('Grouped by number of states, duration >5s ($N_1=%d$, $N_2=%d$ $N_3=%d$ $N_4=%d$)'%tuple(nrsTraces), fontsize=12)
    #ax.set_title('$N_1=%d$, $N_2=%d$ $N_3=%d$ $N_4=%d$'%tuple(nrsTraces), fontsize=12)
    ax.legend(['N=%d'%n for n in nrsTraces[::-1]], fontsize=fontSize)

    return 4



def get_spread(d, center, area):
    
    if np.sum(d)==0:
        return 0,0

    d = 1.0*d/np.sum(d)

    leftOffset = 1
    rightOffset = 1

    nrBins = len(d)

    for m in range(nrBins):
        if (center-leftOffset) >= 0:
            leftBin = d[center-leftOffset]
        else:
            leftBin = 0

        if (center+rightOffset) < nrBins:
            rightBin = d[center+rightOffset]
        else:
            rightBin = 0

        if leftBin > rightBin:
            leftOffset += 1
        elif leftBin < rightBin:
            rightOffset += 1
        else:
            if leftBin==0:
                leftOffset += 1
                rightOffset += 1
            else:
                leftOffset += 1

        cArea = np.sum(d[center-leftOffset:center+rightOffset+1])
        if cArea >= area:
            break

    if (center-leftOffset) < 0:
        leftOffset == center

    if (center+rightOffset) > (nrBins+1):
        rightOffset = nrBins+1-center

    return leftOffset, rightOffset


def get_fret_specs_data(data, descr):
    get_nr_mols = lambda s:int(np.round(np.sum(data[:,descr.index(s)])))

    #nrSelectedMolecules = np.sum([get_nr_mols(i) for i in [
    #    'cat_no_bleach_norm', 
    #    'cat_donor_bleach_only_norm',
    #    'cat_acceptor_bleach_only_norm',
    #    'cat_both_bleached_norm']])
    nrSelectedMolecules = np.sum([get_nr_mols(i) for i in [
        'cat_acceptor_bleach_only_norm',
        'cat_both_bleached_norm']])

    nrZeroMolecules = get_nr_mols('zero_fret')
    nrStableMolecules = get_nr_mols('one_state')
    nrFluctMolecules = np.sum([get_nr_mols(i) for i in 
        ['two_states', 'three_states', 'four_states']])
    
    bins = data[:,descr.index('bins')]
    fret = data[:,descr.index('fret_corr')]
    maxIndex = np.argmax(fret)
    fretMode = bins[maxIndex]
        
    leftWidth, rightWidth = get_spread(fret, maxIndex, 0.75)
    
    leftIndex = maxIndex-leftWidth
    rightIndex = maxIndex+rightWidth
    if leftIndex < 0:
        leftIndex = 0
    if rightIndex >= len(bins):
        rightIndex = len(bins)-1
    fretRange = [fretMode-bins[leftIndex], 
            bins[rightIndex]-fretMode]

    try:
        fractionZero = 1.0*nrZeroMolecules/(nrZeroMolecules+nrSelectedMolecules)
    except:
        fractionZero = 0
    
    try:
        #print nrFluctMolecules, nrStableMolecules
        fractionFluct = 1.0*nrFluctMolecules/(nrStableMolecules+nrFluctMolecules)
    except:
        fractionFluct = 0

    return nrSelectedMolecules, [fretMode, fretRange], [nrZeroMolecules, fractionZero], fractionFluct


def get_signal_hist(data, states, crit, descr='signal'):
    if np.sum(crit) > 0:
        dat = data.select(crit)
        stt = states.select(crit)
        signal = get_state_data(dat, stt, descr=descr)
        cnt, bins, _ = get_fret_histograms(signal, norm_indiv=True, norm_total=False)
        return cnt, bins
    else:
        return np.zeros(HISTOGRAM_NR_BINS), np.linspace(HISTOGRAM_RANGE[0], HISTOGRAM_RANGE[1], HISTOGRAM_NR_BINS)


def plot_stacked_bars(ax, bins, cnt, stacked=True, colors=('b'), *args, **kwargs):
    width = np.mean(np.diff(bins))
    if len(cnt.shape) == 1:
        cnt = np.reshape(cnt, [1,cnt.size])

    nrHist = cnt.shape[0]
 
    p = []
    if not stacked:
        p.append(ax.bar(bins, np.sum(cnt, axis=0), width, align='center', color='w'))
    
    bottom = 0
    for i in range(nrHist):
        if stacked:
            bottom=np.sum(cnt[:i,:],axis=0)
        p.append(ax.bar(bins, cnt[i,:], 
                width=width, 
                bottom=bottom, 
                align='center',
                color=colors[i%len(colors)],
                *args,
                **kwargs))
    return p



def make_allpoints_zero_figure(subsets, panels, ratios):
    #fig = plt.figure(205, figsize=(9,5.5))
    fig = plt.figure(205, figsize=(9,6.5))

    nrPanels = len(subsets)
    nrCols = int(np.ceil(np.sqrt(nrPanels)))
    nrRows = nrCols

    nrFilledRows = np.ceil(1.0*nrPanels/nrCols)
    #fontSize = np.round(6+6.0/nrPanels)
    fontSize=10

    for n in range(len(subsets)):
        data, header = subsets[n]

        try:
            if data==None:
                continue
        except:
            pass

        bins = data[:,0]
        cnt = data[:,9:15]
        descr = header[9:15]
 
        ax = plt.subplot(nrRows, nrCols, n+1)
        
        cnt[np.isnan(cnt)] = 0
        fretCorr = cnt[:,1]
        zeroFret = cnt[:,5]

        plot_stacked_bars(ax, bins, 
                np.vstack((norm_hist(fretCorr)*(1-ratios[n]), norm_hist(zeroFret)*ratios[n])), stacked=True, colors=('w', [0.8,0.8,0.8]))

        ax.set_xlabel('FRET ratio')
        ax.set_ylabel('Count')
        ax.set_xlim(HISTOGRAM_RANGE)
        #ax.legend(['N=%d'%np.sum(fretAll), 'N=%d'%np.sum(zeroFret)], fontsize=fontSize)


        c = n%nrCols
        r = np.floor(n/nrCols)

        ax.set_ylabel('')
        ax.set_yticklabels([])
        ax.set_xticklabels(ax.get_xticks(), fontsize=fontSize+2)
        ax.set_xlabel('FRET ratio', fontsize=fontSize+2)

        crit1 = r < (nrFilledRows-1)
        crit2 = (c >= nrPanels%nrCols) and (r == (nrFilledRows-2)) and (not (nrPanels==(nrFilledRows*nrCols)))
        if crit1:
            if not crit2:
                ax.set_xticklabels([])
                ax.set_xlabel('')
    
    
    # Equalize y-axes
    #for n in range(len(subsets)):
    ymax = max([a.get_ylim()[1] for a in fig.axes])
    [a.set_ylim([0, ymax]) for a in fig.axes]

    return fig
