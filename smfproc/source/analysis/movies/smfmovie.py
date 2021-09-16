""" SMFMOVIE
"""

# System modules
import os, shutil
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import imread as scpimread # deprecated
from imageio import imread as scpimread

# Program modules
from smfproc.source.io.filetypes import MovieFile, TraceFile, SpecsFile
from smfproc.source.io.datatypes import SmfSpecs, SmfSelect
from smfproc.source.io.misc import load_file_paths, search_directory_tree,\
        load_calibration_file, load_laser_profile
from smfproc.source.analysis.movies.environment import LewisLab
from smfproc.source.analysis.movies.analysis import analyze_movie, \
        get_image_background, scale_fluor_intensity


def analyze(args):
    """ Top level movie analysis.
    
    Handles file I/O, data prep for analysis.
    """
    
    filePath, coreSettings, movieSettings = args
    
    mSet = movieSettings
    cSet = coreSettings

    #moviePath = os.path.split(filePath)[0]
    warpMatrix = load_calibration_file(filePath, mSet['REGISTRATION_FILE'])
    if warpMatrix == None:
        return
    
    # Load movie data
    mFile = MovieFile(filePath)
    movie = mFile.read()
    movs = movie.split()

    #movs,_ = load_movie_data(filePath)

    env = LewisLab(mFile)
    #lasers = env.get_all_lasers()
    devices, deviceDescr = env.get_all_devices()
    
    traces, crds = analyze_movie(movs, warpMatrix, mSet)
    nrTraces = traces[0].get_nr_traces()

    laserProfile = load_laser_profile(filePath, mSet['BACKGROUND_FILE'])
    
    #if not (laserProfile==None):
    scale_fluor_intensity(traces[0], traces[1], crds[1], laserProfile)
    
    specs = np.transpose(np.array([
        crds[0][:,0], crds[0][:,1], crds[1][:,0], crds[1][:,1],
        [mSet['MOLECULE_RADIUS']]*nrTraces]))
    specs = SmfSpecs(specs, 
            ['don_row', 'don_col', 'acc_row', 'acc_col', 'radius'])
    
    


    ### === WRITE DATA TO FILE ============================================ ###

    # Create folder for analysis output
    #analysisDir = os.path.splitext(mFile.filePath)[0]
    analysisDir = os.path.splitext(filePath)[0]

    if not os.path.isdir(analysisDir):
        os.mkdir(analysisDir)

    # Write analysis results to file
    header = {'blocks': ['donor', 'acceptor', 'dy', 'dx']}
    traceFilePath = os.path.join(analysisDir, cSet['FILE_BASE']+cSet['TRACE_EXT'])
    traceFile = TraceFile(traceFilePath)
    traceFile.write(traces, header['blocks'])
    
    snapFilePath = os.path.join(analysisDir, 
            cSet['SNAP_BASE']+'1'+cSet['TRACE_EXT'])
    shutil.copyfile(traceFilePath, snapFilePath)

    # Write metadata to file
    devices = [devices.get(n) for n in range(devices.get_nr_traces())]
    metaHeader = {'blocks': deviceDescr}
    traceFile = TraceFile(os.path.join(analysisDir, 
        cSet['FILE_BASE']+cSet['META_EXT']))
    traceFile.write(devices, metaHeader['blocks'])

    # Write coordinates to file
    specsFileName = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['SPECS_EXT'])
    specsFile = SpecsFile(specsFileName)
    specsFile.write(specs)
    #np.savetxt(specsFileName, specs, header='index don_row don_col acc_row acc_col', fmt='%.3f')

    # Clear existing selections file
    selectFileName = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['SELECT_EXT'])
    #if os.path.isfile(selectFileName):
    #    os.remove(selectFileName)

    selectOutFile = SpecsFile(selectFileName)
    #selectOutFile.clear()
    sel = SmfSelect()
    #sel.prefill(['all'], nrTraces)
    sel.set(np.ones(nrTraces), 'all')
    sel.set(np.zeros(nrTraces), 'manual')
    sel.set(np.zeros(nrTraces), 'auto')
    selectOutFile.write(sel)

