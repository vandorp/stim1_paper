""" SMFOUTPUT

Routines for generating analysis output files and figures.
"""

# System packages
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Program modules
from smfproc.source.io.filetypes import SpecsFile, TraceFile
from smfproc.source.io.datatypes import SmfSpecs, SmfSelect, State
from smfproc.source.io.misc import get_output_filepaths as get_filepaths, \
        load_smf_data, save_combined_smf_data

# Program modules - local relative imports
from .analysis import *



def load_hist_data(filePath):
    if not os.path.isfile(filePath):
        return [None, []]

    f = open(filePath, 'r')
    header = f.readline()
    f.close()
    data = np.loadtxt(filePath, float)
    header = [h.strip() for h in header.strip('#').split()]
    return [data, header]

titleFontSize=22
def combine_smf_data_files(dataDirs, outDir):
    trcData = [load_smf_data(d, 'selected') for d in dataDirs]
    save_combined_smf_data(outDir, 'selected', trcData)
    make_analysis_output(outDir, 'selected')
    histData = [load_hist_data(d+'/fret.hist') for d in dataDirs]

    if len(trcData) > 1:
        figs = make_analysis_overview(histData)

        plt.figure(figs[0].number).suptitle(os.path.basename(outDir), fontsize=titleFontSize)
        plt.figure(figs[1].number).suptitle(os.path.basename(outDir), fontsize=titleFontSize)
        plt.figure(figs[2].number).suptitle(os.path.basename(outDir), fontsize=titleFontSize)


        subsetNames = [d.split(os.path.sep)[-1] for d in dataDirs]
        for fig in figs:
            for m in range(len(fig.axes)):
                fig.axes[m].set_title(subsetNames[m])
            plt.figure(fig.number)
            plt.subplots_adjust(top=0.85)

        save_analysis_figure(figs[0].number, 
                os.path.join(outDir, 'overview_allpoints.png'))
        save_analysis_figure(figs[1].number, 
                os.path.join(outDir, 'overview_norm.png'))
        #save_analysis_figure(figs[2].number, 
        #        os.path.join(outDir, 'overview_all_zero.png'))
        save_analysis_figure(figs[2].number, 
                os.path.join(outDir, 'overview_fluct_states.png'))

    titles = [os.path.basename(d) for d in dataDirs]
    fig = make_specs_overview(histData, titles)
    #fig.canvas.set_window_title(os.path.basename(outDir))
    #fig.axes[0].set_title(os.path.basename(outDir))
    save_analysis_figure(fig.number, 
            os.path.join(outDir, 'overview_specs.png'))
    
    #spcs = np.array([SpecsFile(os.path.join(d, 'movie.specs')).read().data
    #        for d in dataDirs])
    spcs = [np.squeeze(SpecsFile(os.path.join(d, 'movie.specs')).read().data)
            for d in dataDirs]
    
    if len(spcs[0])==0:
        return


    spcs = np.array(spcs)
   
    spcs = np.hstack(np.squeeze(spcs))
    spcs = np.reshape(spcs, [spcs.size, 1])
    spcs = SmfSpecs(spcs, ['nr_molecules'])
    sFile = SpecsFile(os.path.join(outDir, 'movie.specs'))
    sFile.write(spcs)


def get_smf_data_subset(data, selectionId):
    traces = data[0]
    if traces==None:
        return

    traces, header = data[0]
    if len(header) == 0:
        return
    
    specs = data[1]
    sel = data[2]

    crit = SmfSelect(sel).get(selectionId)

    if crit.size == 0:
        return [], None, None, None
            
    #traces, trcHeader = trc
    traces = [trc.select(crit) for trc in traces]
    specs = specs.select(crit)
    sel = sel.select(crit)

    data = [[traces, header], specs, SmfSelect(sel), None]
    return data


def save_movie_specs(baseDir, fileName):
    sFile = SpecsFile(os.path.join(baseDir, fileName))
    specs = sFile.read()
    
    nrDetected = specs.data.shape[0]

    mSpecs = SmfSpecs(nrDetected, ['nr_molecules'])
    mFile = SpecsFile(os.path.join(baseDir, 'movie.specs'))
    mFile.write(mSpecs)


def analyze_movie_traces(args):
    movieFilePath, coreSettings, outputSettings, manual = args
    basePath = os.path.splitext(movieFilePath)[0]

    movieDescr = os.path.sep.join(movieFilePath.split(os.path.sep)[-3:])
    print('Producing analysis outputs for movie %s'%movieDescr)

    sourceFlag = 'detected'

    # Export movie specs summary
    save_movie_specs(basePath, sourceFlag+'.specs')
    
    # Export traces subset
    data = load_smf_data(basePath, sourceFlag) 
    
    if manual:
        selDescr = ['manual']
    else:
        selDescr = ['selected']
    
    for descr in selDescr:
        subset = get_smf_data_subset(data, descr)
        save_combined_smf_data(basePath, 'selected', [subset])
    
    # Make analysis figures
    make_analysis_output(basePath, sourceFlag, manual)


def make_analysis_output(basePath, flag, manual=False):

    data = load_smf_data(basePath, flag)

    ftitle = os.path.basename(basePath)
    
    # === Traces === #
    traces = data[0]
    if traces==None:
        return

    trc, header = data[0]
    if len(header) == 0:
        return


    # === Specs === #
    gamma = data[1].get('gamma')
    selections = data[2]

    # === Output === #
    fret = trc[header.index('fret')]
    states = State(trc[header.index('states')])


    if manual:
        selection_type = 'manual'
    else:
        selection_type = 'auto'
    
    cnt1, bins, descr1 = get_categories_histograms_data(fret, states, selections, selection_type)
    cnt2, bins, descr2 = get_fret_histograms_data(fret, states, selections, gamma, selection_type)
    cnt3, bins, descr3 = get_nrstates_histograms_data(fret, states, selections, gamma, selection_type)
    cnt4, bins, descr4 = get_hmm_histograms_data(fret, states, selections, gamma, selection_type)
    
    cnt = np.vstack((cnt1, cnt2, cnt3, cnt4))
    data = np.transpose(np.vstack((bins, cnt)))
    histHeader = 'bins '+' '.join(descr1+descr2+descr3+descr4)

    np.savetxt(os.path.join(basePath, 'fret.hist'), data, '%.3f', header=histHeader)

    # === Output for ebFRET === #
    donor = trc[header.index('donor')]
    acceptor = trc[header.index('acceptor')]
    stackedData = stack_fret_traces(donor, acceptor, states, selections, gamma)
    np.savetxt(os.path.join(basePath, 'fret.dat'), stackedData, fmt='%.7e', delimiter='\t')


    fig1 = plt.figure(100, figsize=(6,9))
    fig2 = plt.figure(101)
    ax2 = plt.gca()
    fig3 = plt.figure(102)
    ax3 = plt.gca()
    fig4 = plt.figure(103)
    ax4 = plt.gca()
    fig5 = plt.figure(104)
    ax5 = plt.gca()
    #fig6 = plt.figure(105)
    #ax6 = plt.gca()
    #fig7 = plt.figure(106)
    #ax7 = plt.gca()

    #ax2.set_title(ftitle)
    #ax3.set_title(ftitle)
    #ax4.set_title(ftitle)
    #ax5.set_title(ftitle)

    fig1.suptitle('fret categories ('+selection_type+')', fontsize=titleFontSize)
    fig2.suptitle('fret all points ('+selection_type+')', fontsize=titleFontSize)
    fig3.suptitle('fret indiv. norm ('+selection_type+')', fontsize=titleFontSize)
    fig4.suptitle('fret nr. states ('+selection_type+')', fontsize=titleFontSize)
    fig5.suptitle('fret hmm model ('+selection_type+')', fontsize=titleFontSize)
    #fig6.suptitle(ftitle, fontsize=titleFontSize)
    #fig7.suptitle(ftitle, fontsize=titleFontSize)


    get_categories_histograms_figure(fig1, cnt1, bins, descr1)
    get_fret_histograms_figures(ax2, ax3, np.transpose(cnt2), bins, descr2)
    get_nrstates_histograms_figure(ax4, np.transpose(cnt3), bins, descr3)
    get_hmm_histograms_figures(ax5, np.transpose(cnt4), bins, descr4)
    #get_fret_histograms_figures_manual(ax6, ax7, np.transpose(cnt2), bins, descr2)

    save_analysis_figure(100, os.path.join(basePath, 'fret_categories.png'))
    save_analysis_figure(101, os.path.join(basePath, 'fret_allpoints.png'))
    save_analysis_figure(102, os.path.join(basePath, 'fret_norm.png'))
    save_analysis_figure(103, os.path.join(basePath, 'fret_fluct_states.png'))
    save_analysis_figure(104, os.path.join(basePath, 'fret_model.png'))
    #save_analysis_figure(105, os.path.join(basePath, 'fret_allpoints_manual.png'))
    #save_analysis_figure(106, os.path.join(basePath, 'fret_norm_manual.png'))




def make_analysis_overview(subsets):
    
    fig1 = plt.figure(200)
    fig2 = plt.figure(201)
    fig3 = plt.figure(202)
    fig4 = plt.figure(203)

    nrPanels = len(subsets)
    nrCols = int(np.ceil(np.sqrt(nrPanels)))
    #nrRows = (nrPanels%nrCols)+1
    nrRows = nrCols

    nrFilledRows = np.ceil(1.0*nrPanels/nrCols)

    fontSize = np.round(6+6.0/nrPanels)

    for n in range(len(subsets)):
        data, header = subsets[n]

        try:
            if data==None:
                continue
        except:
            if len(data)==0:
                continue

        bins = data[:,0]
        cnt1 = data[:,9:15]
        descr1 = header[9:15]
        cnt2 = data[:,15:19]
        descr2 = header[15:19]
 
        fig1 = plt.figure(200)
        ax1 = plt.subplot(nrRows, nrCols, n+1)
        fig2 = plt.figure(201)
        ax2 = plt.subplot(nrRows, nrCols, n+1)
        fig3 = plt.figure(202)
        ax3 = plt.subplot(nrRows, nrCols, n+1)
        #fig4 = plt.figure(203)
        #ax4 = plt.subplot(nrRows, nrCols, n+1)
       
        get_fret_histograms_figures(ax1, ax2, cnt1, bins, descr1, fontSize=fontSize)
        get_nrstates_histograms_figure(ax3, cnt2, bins, descr2, fontSize=fontSize)

        c = n%nrCols
        r = np.floor(n/nrCols)
        
        for ax in [ax1, ax2, ax3]:
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
    for fig in [fig1, fig2, fig3]:
        ymax = max([a.get_ylim()[1] for a in fig.axes])
        [a.set_ylim([0, ymax]) for a in fig.axes]

    return [fig1, fig2, fig3]



def spread_plot(ax, x, y, prange, width=0.8, *args, **kwargs):
    nrBars = len(x)

    for n in range(nrBars):

        origin = (x[n]-width/2, y[n]-prange[0][n])
        height = prange[0][n]+prange[1][n]
        ax.add_patch(
                patches.Rectangle(
                    origin,
                    width,
                    height,
                    *args,
                    **kwargs))
        ax.plot([origin[0], origin[0]+width], [y[n], y[n]], 'k')

    
def get_all_specs(subsets):
    allSelected = []
    allZero = []
    allFluct = []
    allFret = []
   
    for k in range(len(subsets)):
        data, header = subsets[k]

        try:
            if data==None:
                selected = fluct = 0
                zero = fret = [np.nan, [0,0]]
            else:
                selected, fret, zero, fluct = get_fret_specs_data(data, header)
        except:
            selected, fret, zero, fluct = get_fret_specs_data(data, header)

        allSelected.append(selected)
        allZero.append(zero)
        allFret.append(fret)
        allFluct.append(fluct)
    
    return allSelected, allZero, allFret, allFluct


def make_specs_overview(subsets, titles):
    
    allSelected, allZero, allFret, allFluct = get_all_specs(subsets)
    
    #fig1 = plt.figure(200)
    
    #nrPanels = len(subsets)
    #nrCols = int(np.ceil(np.sqrt(nrPanels)))
    #nrRows = (nrPanels%nrCols)+1
    #nrRows = nrCols

    #nrFilledRows = np.ceil(1.0*nrPanels/nrCols)

    #fontSize = np.round(6+6.0/nrPanels)
    
    nrSubsets = len(subsets)
    
    if nrSubsets > 2:
        figWidth = nrSubsets
    else:
        figWidth = 3

    fig = plt.figure(figsize=(figWidth,9))
    #fig = plt.figure()
    ax = plt.subplot(411)
    ax.bar(np.arange(len(subsets))+0.6, allSelected, width=0.8, color='r')
    ax.set_xlim([0, len(subsets)+1])
    ax.set_ylabel('count')
    ax.set_title('Selected molecules')
    ax.set_xticks(range(1,len(subsets)+1))
    ax.set_xticklabels([])

    ax = plt.subplot(412)
    ax.plot([0,nrSubsets+1], [0,0], '--')
    spread_plot(ax, 
            np.arange(len(subsets))+1, 
            [f[0] for f in allFret], 
            [[f[1][0] for f in allFret], [f[1][1] for f in allFret]],
            color=[0.7,0.7,0.7])
    #ax.plot(np.arange(len(subsets))+1, [f[0] for f in allFret], 'o', color=[0.7,0.7,0.7])
    #ax.errorbar(
    #        np.arange(len(subsets))+1, 
    #        [f[0] for f in allFret], 
    #        yerr=[[f[1][0] for f in allFret], [f[1][1] for f in allFret]],
    #        fmt='o', color=[0.7,0.7,0.7])

    ax.set_xlim([0, len(subsets)+1])
    ax.set_ylim([-0.25,1.25])
    ax.set_ylabel('mode')
    ax.set_title('FRET ratio')
    ax.set_xticks(range(1,len(subsets)+1))
    ax.set_xticklabels([])
    
    ax = plt.subplot(413)
    #ax.bar(np.arange(len(subsets))+0.6, np.array([z[0] for z in allZero]), width=0.8, color='g', alpha=0.25)
    #ax.bar(np.arange(len(subsets))+0.6, np.array([z[1] for z in allZero]), width=0.8, color='g')
    ax.set_xlim([0, len(subsets)+1])
    ax.set_ylim([0,1.1])
    ax.set_ylabel('fraction')
    ax.set_title('Zero-FRET molecules')
    ax.set_xticks(range(1,len(subsets)+1))
    ax.set_xticklabels([])

    ax = plt.subplot(414)
    ax.bar(np.arange(len(subsets))+0.6, allFluct, width=0.8, color='b')
    ax.set_xlim([0, len(subsets)+1])
    ax.set_ylim([0,1.1])
    ax.set_ylabel('fraction')
    ax.set_title('Fluctuating molecules')
    ax.set_xticks(range(1,len(subsets)+1))
    ax.set_xticklabels(titles, rotation=45, ha='right', fontsize=8)
    
    plt.tight_layout()

    return fig 
