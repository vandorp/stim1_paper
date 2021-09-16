""" SMFTRACES
"""

# System modules
import os, shutil
import numpy as np

# Program modules - relative imports
from .analysis import get_signal_states, correct_xtalk, correct_baseline, \
        get_fret_ratio, get_initial_trace_selection

# Program modules
from smfproc.source.io.filetypes import TraceFile, SpecsFile
from smfproc.source.io.datatypes import SmfSpecs, SmfSelect
from smfproc.source.io.misc import load_file_paths, collect_trace_data


def analyze(args):
    """ Top level raw trace analysis.
    
    Handles file I/O, some data prep.
    """

    movieFilePath, coreSettings, traceSettings = args

    tSet = traceSettings
    cSet = coreSettings

    analysisDir = os.path.splitext(movieFilePath)[0]

    sourceFilePath = os.path.join(analysisDir, 
            cSet['SNAP_BASE']+'1'+cSet['TRACE_EXT'])
    metaFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['META_EXT'])
    selectFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['SELECT_EXT'])
    outputFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['TRACE_EXT'])
    snapFilePath = os.path.join(analysisDir, 
            cSet['SNAP_BASE']+'2'+cSet['TRACE_EXT'])

    sourceFile = TraceFile(sourceFilePath)
    metaFile = TraceFile(metaFilePath)
    selectFile = SpecsFile(selectFilePath)
    outputFile = TraceFile(outputFilePath)

    traces, tracesHeader = sourceFile.read()
    lasers, lasersHeader = metaFile.read()

    traces = traces + lasers
    descr = tracesHeader + lasersHeader

    don = traces[descr.index('donor')]
    acc = traces[descr.index('acceptor')]
    
    # Determine states and make corrections
    correct_xtalk(don, acc, tSet['XTALK'])

    states = get_signal_states(traces, descr,
            tSet['FILTER_WINDOW'],
            tSet['SAVGOL_ORDER'],
            tSet['SIGNAL_SNR'],
            tSet['BLEACH_WINDOW'],
            tSet['SIGNAL_WINDOW'])

    correct_baseline(don, acc, states)

    fret = get_fret_ratio(don, acc)
    
    # Write traces to file    
    outputFile.write_block(don, 'donor')
    outputFile.write_block(acc, 'acceptor') 
    outputFile.write_block(fret, 'fret')
    outputFile.write_block(states, 'states')
    
    shutil.copyfile(outputFilePath, snapFilePath)

    # === Write selections to file ===
    selections = get_initial_trace_selection(don, acc, states,
            tSet['FILTER_WINDOW'],
            tSet['STATE_SNR'])

    # The selections file currently has only the 'all' and 'manual' columns
    currentSelections = SmfSelect(selectFile.read())
    
    # Add the new selections
    currentSelections.merge(selections)

    # Make final selections
    #autoRejected = currentSelections.get(['no_fret', 'no_initial_signal'], np.logical_or)
    autoRejected = currentSelections.get('no_fret')
    currentSelections.set(autoRejected==False, 'auto')

    finalSelected = currentSelections.get(['manual', 'auto'], np.logical_xor)
    currentSelections.set(finalSelected, 'selected')
    
    # Write to file
    selectFile.write(currentSelections)
