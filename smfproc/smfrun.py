""" SMFRUN.py - Batch processing of single-molecule FRET data
=============================================================

* Handles uManager stacked-tiff movies (.tif) and a legacy custom format (.pma)
* Assumes movies have a vertical split with donor channel on the left
* Requires a donor-acceptor channel registration file named 'calibration.txt'
to be present in the data folder

Usage:

Issuing 'python smfrun.py' starts the GUI.

'python smfrun.py --path [folder]'
Starts the GUI and loads the movies in the specified folder.

'python smfrun.py --path [folder] --analyze'
Full analysis of the movies in specified folder, without opening the GUI.

'python smfrun.py --path [folder] --analyze [step]'
Partial analysis of the movies in specified folder, starting from [step].
[step] can be any of the five options above. E.g. if you want to change some trace
selection parameters, but don't want to redo all the trace fitting (the most time-
consuming part) use --analyze states.

'python smfrun.py --help'
Print all command line options.

For use with Python 3. Cross-platform in principle. Developed in Ubuntu 16 using 
Python 3.5. Somewhat tested on Windows 10, basically untested on Mac.

Supplement to
"Conformational dynamics of auto-inhibition in the ER calcium sensor STIM1"
Stijn van Dorp, Ruoyi Qiu, Ucheor B. Choi, Minnie M. Wu, Michelle Yen, 
Michael Kirmiz, Axel T. Brunger, Richard S. Lewis
https://doi.org/10.1101/2020.12.17.423361
"""


""" ===========================================================================
= IMPORTS
=========================================================================== """

# System modules
import os, shutil, glob
import numpy as np
from argparse import ArgumentParser

# Program modules
from smfproc.source.io.misc import load_file_paths
from smfproc.source.io.filetypes import AnalysisFile
from smfproc.source.batch.smfbatch import run_batch_analysis
from smfproc.source.gui.gui_setup import build_gui


""" ===========================================================================
= USER-DEFINED SETTINGS
=========================================================================== """

# Average the first x frames of the movie for signal detection.
MOVIE_SIGNAL_FRAMES = 5

# Threshold for signal detection. Multiplied by the std of the pixel values
# of the averaged signal frames.
MOVIE_PIXEL_SNR = 1.5

# For signal integration. Diameter is 2xradius + 1
MOVIE_MOLECULE_RADIUS = 2

# Area radius for local dynamic background subtraction.
MOVIE_BACKGROUND_RADIUS = 17

# Donor to acceptor channel cross-talk.
TRACE_XTALK = 0.07

# Maximum number of states allowed after HMM fitting of the FRET signal. Traces
# with more than this number of distinct states will not be included.
TRACE_MAX_NR_STATES = 4

# Criteria for trace selection
SELECT_REQUIRE_ACCEPTOR_BLEACH = True
SELECT_REQUIRE_DONOR_BLEACH = False
SELECT_GAMMA_RANGE = [0.5, 2.5]
SELECT_MINIMUM_SNR = 5




""" ===========================================================================
*******************************************************************************
*** BELOW THIS LINE NOTHING SHOULD NOT NEED TO BE ALTERED *********************
*******************************************************************************
=========================================================================== """




""" ===========================================================================
= DEFINITIONS
=========================================================================== """

def wrap_settings(varList, flag=''):
    """ Collects and evaluates named variables.

    Variables have to exist in the current workspace. Variable names containing
    the provided flag are evaluated and stored in a dictionary.
    
    Arguments
    ---------
        varList : list of variable names
        flag : string 

    Returns
    -------
        Dictionary with evaluated variables
     """

    hasFlag = np.array([v.find(flag)==0 for v in varList])
    names = [v for v in np.array(varList)[hasFlag]]
    shortNames = [v.split(flag)[1].lstrip('_') for v in np.array(varList)[hasFlag]]

    settings = dict()
    for n in range(len(names)):
        settings[shortNames[n]] = eval(names[n])
    
    return settings


def cmd_args():
    """ Parses command line arguments.
    
    Returns
    -------
        Dictionary with values of command line arguments.

    """

    usage = "%(prog)s [options]"
    parser = ArgumentParser(
            description='View, analyze and export smf movies, traces, and model fits.', 
            usage=usage)
    parser.add_argument('-p', '--path', type=str,
            action='store',
            default=None,
            help="Directory or file to analyze or view.")
    parser.add_argument('-s', '--dataset', type=str,
            action='store',
            default=None,
            help="Dataset to analyze. Need to also specify the data root folder using --path option.")
    parser.add_argument('-a', '--analyze', 
            action='store',
            nargs='?',
            const='laser',
            choices=('laser', 'movies', 'traces', 'states', 'output'),
            help="Run analysis")
    parser.add_argument('-x', '--purge',
            action='store_true',
            help="Delete all existing analysis files")
    parser.add_argument('-u', '--update',
            action='store_true',
            help="Only analyze if no previous analysis files exist")
    parser.add_argument('--nomovie',
            action='store_true',
            help="Run GUI without loading movies")
    parser.add_argument('-w', '--stepwise',
            action='store_true',
            help="Only run selected analysis step, skip subsequent analyses")
    parser.add_argument('-d', '--nrcores', type=int,
            action='store',
            default=0,
            help="Number of cores to use for parallel analyses. Default behavior is to use all available cores, except for one.")
    parser.add_argument('-m', '--manual',
            action='store_true',
            help="Only export manually selected traces for analysis.")

    args = parser.parse_args()
    return vars(args)



""" ===========================================================================
= ANALYSIS SETTINGS
=========================================================================== """

# Flags and naming conventions for input and output files
MOVIE_EXTENSIONS = ['.pma', '.tif']
MOVIE_REGISTRATION_FILE = 'calibration.txt'
MOVIE_BACKGROUND_FILE = 'laser_profile.png'

CORE_FILE_BASE = 'detected'
CORE_SNAP_BASE = 'snap'
CORE_TRACE_EXT = '.trc'
CORE_META_EXT = '.meta'
CORE_SPECS_EXT = '.specs'
CORE_SELECT_EXT = '.select'

OUTPUT_SUBSET = 'selected'

# Movie analysis settings
MOVIE_DONOR_LEFT = True
MOVIE_BACKGR_RESOLUTION = 28

# Savitzky-Golay filter. This is a smoothing average (convolution)
# filter that preserves features such as sharp transitions. Used for
# detecting bleaching events.
TRACE_FILTER_WINDOW = 9
TRACE_SAVGOL_ORDER = 3

TRACE_BLEACH_WINDOW = 3
TRACE_SIGNAL_WINDOW = 5

# Thresholds for signal detection. Multiplied by the std of the noise of the 
#filtered trace.
TRACE_SIGNAL_SNR = 3.0
TRACE_STATE_SNR = 2.0
STATE_TRANSITION_SNR = 2.0
STATE_BLEACH_SNR = 0.25


""" ===========================================================================
= MAIN PROGRAM
=========================================================================== """

# Collect variables in the workspace that contain program settings.
allVars = dir()

settings = dict()
settings['core'] = wrap_settings(allVars, 'CORE')
settings['movie'] = wrap_settings(allVars, 'MOVIE')
settings['trace'] = wrap_settings(allVars, 'TRACE')
settings['state'] = wrap_settings(allVars, 'STATE')
settings['select'] = wrap_settings(allVars, 'SELECT')
settings['output'] = wrap_settings(allVars, 'OUTPUT')


def main():
    """ Main program.

    Based on command line arguments this either runs analysis or starts the GUI.
    If the --analyze option was set, analysis will run in the command line. If it
    was omitted, the GUI will be started.

    Analysis can be directed to a single data folder indicated by the --path 
    argument, or to groups of folders specified in a dataset file, and indicated
    by the --dataset argument.

    Returns
    -------
        
    """
    
    # Parse command line arguments. Run --help for descriptions of the various
    # options.
    cmdArgs = cmd_args()
    
    argPath = cmdArgs['path']
    dataset = cmdArgs['dataset']
    runAnalysis = cmdArgs['analyze']
    purge = cmdArgs['purge']
    update = cmdArgs['update']
    loadMovie = cmdArgs['nomovie'] == False
    cascade = cmdArgs['stepwise'] == False
    distributed = cmdArgs['nrcores']
    manual = cmdArgs['manual']

    launchGUI = (not runAnalysis) and (not dataset)
    

    if dataset:
        # Executed if a path to a dataset file was provided by
        # the --dataset option.
        #
        # Dataset files contain groups of paths to data folders
        # that belong together, i.e. were derived from the same
        # sample. The --path argument needs to contain the root
        # path if the dataset file contains relative file paths.

        anFile = AnalysisFile(dataset)
        groups, dataDirs, header = anFile.read()
        allDataSubs = []
        [allDataSubs.extend(d) for d in dataDirs]
        allDataDirs = []
        [allDataDirs.extend(d) for d in allDataSubs]
        allDataDirs = [os.path.join(argPath, d) for d in allDataDirs]
    else:
        # If no dataset file was provided, the --path argument
        # points directly to a data folder.
        allDataDirs = [argPath]

    if runAnalysis:
        # Run analysis if --analyze was set. Runs either a full or partial
        # analysis, depending on the argument of --analyze and the --stepwise
        # setting.

        for dDir in allDataDirs:
            shortDir = os.path.sep.join(dDir.split(os.path.sep)[-3:])
            print('>>> Batch analysis on %s <<<'%shortDir)

            moviePaths = load_file_paths(dDir, None, MOVIE_EXTENSIONS)
            run_batch_analysis(moviePaths, runAnalysis, 
                    purge, update, cascade, distributed, manual, settings)
            
    if launchGUI:
        # The GUI is launched if --analyze was not specified. Data files or 
        # folders provided by --path or --dataset will be loaded into the GUI
        # for inspection and/or analysis.

        gui = build_gui()
        gui.settings = settings

        movieGUI = gui.get_item('frames').get_item('movie')
        if loadMovie:
            movieGUI.get_var('show_movie').set("First frames only")
        else:
            movieGUI.get_var('show_movie').set("Skip movie")

        inputPath = allDataDirs[0]
        if inputPath == None:
            moviePaths = []
        else:
            if os.path.splitext(inputPath)[-1] == '.trc':
                moviePaths = [inputPath]
            else:
                moviePaths = load_file_paths(inputPath, None, MOVIE_EXTENSIONS)

        batchGUI = gui.get_item('movies').get_item('batch_sub')
        batchGUI.update(moviePaths)
        gui.mainloop()



if __name__ == "__main__":
    """ Run program if executed from the command line. """

    main()
