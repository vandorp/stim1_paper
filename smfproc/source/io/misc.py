""" Miscellaneous file searching, reading, and writing. """


# System packages
import os, fnmatch
import numpy as np

# Program modules
from smfproc.source.io.filetypes import TraceFile, SpecsFile
from smfproc.source.io.datatypes import SmfSelect


def load_calibration_file(filePath, calibrationFile):
    if (calibrationFile==None) or (calibrationFile==''):
        print('No image calibration provided. Aborting.')
        return None

    # Get image calibration data
    imageRegisterFile = search_directory_tree(
            filePath, calibrationFile)

    if (imageRegisterFile==None):
        print('No image calibration found. Aborting.')
        return None
    else:
        printOut = os.path.join(*imageRegisterFile.split(os.path.sep)[-3:])
        print('Image calibration found: %s'%printOut)

    warpMatrix = np.vsplit(np.loadtxt(imageRegisterFile), 2)
    return warpMatrix


def load_laser_profile(filePath, profileFileName):

    # Get image calibration data
    profileFilePath = search_directory_tree(
            filePath, profileFileName)

    if (profileFilePath==None):
        print('No laser profile found.')
        return None
    else:
        printOut = os.path.join(*profileFilePath.split(os.path.sep)[-3:])
        print('Laser profile found: %s'%printOut)

    return plt.imread(profileFilePath)


def search_directory_tree(dirPath, flag):
    files = []
    #pathDirs = ['.']+dirPath.split(os.path.sep)[:-1]
    pathDirs = dirPath.split(os.path.sep)[:-1]

    for p in range(len(pathDirs)+1):
        pathDir = os.path.sep.join(pathDirs[:p])
        filePath = pathDir+os.path.sep+'*'+flag+'*'
        files.extend(glob.glob(filePath))

    if len(files) > 0:
        return files[-1]
    else:
        return None




def get_dir_files(dirName, strExpr='*.*', recurs=True):
    """Recursively search directory and subdirectories for filenames with a
    specific pattern."""

    matches = []
    for root, dirnames, filenames in os.walk(dirName):
        for filename in fnmatch.filter(filenames, strExpr):
            matches.append(os.path.abspath(os.path.join(root, filename)))

        if not recurs:
            # The 'dirnames' list holds the remaining subdirectories
            del dirnames[:]
    return sorted(matches)


def load_file_paths(dirPath, fileName, flags):
    if dirPath==None:
        return []

    if fileName==None:
        filePaths = [get_dir_files(dirPath, '*'+flag) for flag in flags]
        filePaths = sum(filePaths, [])
    else:
        filePaths = [os.path.join(dirPath, fileName)]
    return filePaths


def get_output_filepaths(fPath, cSet, snapIndex):
    analysisDir = os.path.splitext(fPath)[0]

    sourceFilePath = os.path.join(analysisDir, 
            cSet['SNAP_BASE']+str(snapIndex)+cSet['TRACE_EXT'])
    selectFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['SELECT_EXT'])
    specsFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['SPECS_EXT'])
    outputFilePath = os.path.join(analysisDir, 
            cSet['FILE_BASE']+cSet['TRACE_EXT'])
    snapFilePath = os.path.join(analysisDir, 
            cSet['SNAP_BASE']+str(snapIndex+1)+cSet['TRACE_EXT'])

    return sourceFilePath, selectFilePath, specsFilePath, outputFilePath, snapFilePath


def collect_trace_data(tracePath, cSet):
    fileBase = os.path.splitext(tracePath)[0]

    metaPath = fileBase+cSet['META_EXT']
    specsPath = fileBase+cSet['SPECS_EXT']
    selectPath = fileBase+cSet['SELECT_EXT']

    traceFile = TraceFile(tracePath)

    traces, tracesHeader = traceFile.read()
    lasers, lasersHeader = TraceFile(metaPath).read()
    specs = SpecsFile(specsPath).read()
    selections = SpecsFile(selectPath).read()
    selections = SmfSelect(selections)

    traces = traces + lasers
    descr = tracesHeader + lasersHeader

    if 'time' in descr:
        timeInd = descr.index('time')
        for ts in traces:
            time = traces[timeInd].data
            ts.time = np.reshape(time, time.size)

    return traces, descr, specs, selections


def load_smf_data(dataDir, flag):
    exts = ['.trc', '.specs', '.select', '.meta']
    filePaths = [os.path.join(dataDir, flag+ext) for ext in exts]
            
    fileTypes = [TraceFile, SpecsFile, SpecsFile, TraceFile]
    
    outData = []
    for n in range(4):
        filePath = filePaths[n]
        fileType = fileTypes[n]
        if os.path.isfile(filePath):
            try:
                outData.append(fileType(filePath).read())
            except:
                outData.append(None)
        else:
            print('*** File %s does not exist ***'%tuple(exts)[n])
            outData.append(None)

    return outData


def save_combined_smf_data(outDir, flag, data):
    tracesOutFilePath = os.path.join(outDir, flag + '.trc')
    tracesOutFile = TraceFile(tracesOutFilePath)
    tracesOutFile.clear()
    
    specsOutFilePath = os.path.join(outDir, flag + '.specs')
    specsOutFile = SpecsFile(specsOutFilePath)
    specsOutFile.clear()
    
    selectOutFilePath = os.path.join(outDir, flag + '.select')
    selectOutFile = SpecsFile(selectOutFilePath)
    selectOutFile.clear()

    for d in data:
        traces, specs, select, meta = d
        if traces==None:
            continue

        trc, header = traces

        if len(header) > 0:
            tracesOutFile.append(trc, header)
            specsOutFile.append(specs)
            selectOutFile.append(select)
            # Skip meta for now



