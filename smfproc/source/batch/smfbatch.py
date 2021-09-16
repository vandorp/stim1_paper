""" SMFBATCH """

# Python packages
import os, shutil
import numpy as np
import multiprocessing as mp
import scipy.misc as scymisc

# Program modules
from smfproc.source.analysis.movies.smfmovie import analyze as analyze_movie
from smfproc.source.analysis.traces.smftraces import analyze as analyze_traces
from smfproc.source.analysis.states.smfstates import analyze as analyze_states, refine as refine_states
from smfproc.source.analysis.output.smfoutput import analyze_movie_traces, combine_smf_data_files
from smfproc.source.io.misc import load_smf_data

# Program modules - local relative imports
from .batch_analysis import generate_laser_profile
from .batch_misc import get_shared_base_path

def run_batch(batchArgs):
    func, args = batchArgs
    [func(arg) for arg in args]


def distribute(distrArgs): 
    func, args, pool = distrArgs
    nrCores = pool._processes
    distrArgs = [[func, args[n::nrCores]] for n in range(nrCores)]
    pool.map(run_batch, distrArgs)


def run_batch_analysis(moviePaths, runAnalysis, purge, update, cascade, distributed, manual,
        settings):
    

    analysisDirs = [os.path.splitext(p)[0] for p in moviePaths]

    analyzeMovie = np.repeat(True, len(moviePaths))
    for n in range(len(analysisDirs)):
        adir = analysisDirs[n]

        if os.path.isdir(adir):
                        
            snapFile = os.path.join(adir, 'snap4.trc')
            if update:
                if os.path.isfile(snapFile):
                    analyzeMovie[n] = False
                else:
                    if purge:
                        shutil.rmtree(adir)
            else:
                if purge:
                    shutil.rmtree(adir)

    
    moviePaths = [moviePaths[i] for i in np.where(analyzeMovie)[0]]
    nrMovies = len(moviePaths)

    if (runAnalysis==None) or (nrMovies==0):
        return
    
    print('*** Analyzing %d movies ***'%nrMovies)

    generateLaserProfile = runAnalysis=='laser'
    analyzeMovies = (runAnalysis=='movies') | (generateLaserProfile & cascade)
    analyzeTraces = (runAnalysis=='traces') | (analyzeMovies & cascade)
    analyzeStates = (runAnalysis=='states') | (analyzeTraces & cascade)
    outputAnalysis = (runAnalysis=='output') | (analyzeStates & cascade)

    nrAvailCores = mp.cpu_count()
    if distributed <= 0:
        nrCores = nrAvailCores-1
    elif distributed > nrAvailCores:
        nrCores = nrAvailCores
    else:
        nrCores = distributed
    
    if nrCores > 1:
        pool = mp.Pool(processes=nrCores)
    
    print('*** Processing on %d cores ***'%nrCores)

    if generateLaserProfile:
        movBckFile = settings['movie']['BACKGROUND_FILE']
        movBckRes = settings['movie']['BACKGR_RESOLUTION']
        
        sharedPath = get_shared_base_path(moviePaths)[:-1]
        argPath = os.path.sep.join(sharedPath)
        outputFile = os.path.join(argPath, movBckFile)
        
        laserProfile = generate_laser_profile(moviePaths, movBckRes)
        try:
            scymisc.imsave(outputFile, laserProfile)
        except:
            raise ImportError('Package IMSAVE not available for generating a laser profile.')


    if analyzeMovies:
        args = [[mPath, settings['core'], settings['movie']] 
                for mPath in moviePaths]

        if nrCores > 1:
            distribute([analyze_movie, args, pool])
        else:
            run_batch([analyze_movie, args])

    if analyzeTraces:
        args = [[mPath, settings['core'], settings['trace']] 
                for mPath in moviePaths]

        if nrCores > 1:
            distribute([analyze_traces, args, pool])
        else:
            run_batch([analyze_traces, args])

        if nrCores > 1:
            distribute([analyze_states, args, pool])
        else:
            run_batch([analyze_states, args])

    if analyzeStates:
        args = [[mPath, settings['core'], settings['trace'], settings['state'], \
                settings['select']] for mPath in moviePaths]
        if nrCores > 1:
            distribute([refine_states, args, pool])
        else:
            run_batch([refine_states, args])

    if outputAnalysis:
        args = [[mPath, settings['core'], settings['output'], manual] 
                for mPath in moviePaths]
        #if distributed:
        #    distribute([analyze_movie_traces, args, pool])
        #else:
        run_batch([analyze_movie_traces, args])
        
        if len(moviePaths)>1:
            inputDirs = [os.path.splitext(m)[0] for m in moviePaths]
            combineDir = os.path.sep.join(get_shared_base_path(moviePaths)[:-1])
            combine_smf_data_files(inputDirs, combineDir)
        
