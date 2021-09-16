""" Locating and exporting files for batch analysis. """


# System packages
import os
import numpy as np

# Program modules
from smfproc.source.io.filetypes import TraceFile, SpecsFile
from smfproc.source.io.datatypes import SmfSelect



build_source = lambda s,d: s['FILE_BASE']+s[d]
build_path = lambda M,D: os.path.join(os.path.splitext(M)[0], D)




def get_shared_base_path(paths):
    sharedPath = []
    if len(paths)==0:
        return sharedPath
    
    #nrChars = len(movieFiles[0])
    C = True
    n = 0
    dirs = [f.split(os.path.sep) for f in paths]
    while C:
        for m in range(len(paths)-1):
            try:
                currDirs = dirs[m]
                nextDirs = dirs[m+1]

                currDir = currDirs[n]
                nextDir = nextDirs[n]

                if not (currDir == nextDir):
                    C = False
                    break
            except:
                C = False
                break
        
        if len(paths) > 1:
            sharedPath.append(currDir)
            n += 1
        else:
            C = False
            sharedPath = paths[0].split(os.path.sep)
            
    return sharedPath


def get_traces_subset(mFile, sources, selectionId):
    fSel = build_path(mFile, build_source(sources, 'SELECT_EXT'))
    fTrc = build_path(mFile, build_source(sources, 'TRACE_EXT'))
    fSpc = build_path(mFile, build_source(sources, 'SPECS_EXT'))
    
    sel = SpecsFile(fSel).read()
    crit = SmfSelect(sel).get(selectionId)

    if crit.size == 0:
        return [], None, None
            
    traces, trcHeader = TraceFile(fTrc).read()
    traces = [trc.select(crit) for trc in traces]

    specs = SpecsFile(fSpc).read()
    specs = specs.select(crit)

    sel = sel.select(crit)

    return [traces, trcHeader], specs, sel



def export_traces_subset(selectId, movieFiles, settings):

    sharedPath = get_shared_base_path(movieFiles)
    basePath = os.path.sep.join(sharedPath[:-1])

    #settings = self.main.settings['core']
    tracesOutFilePath = \
            os.path.join(basePath, selectId+settings['TRACE_EXT'])
    specsOutFilePath = \
            os.path.join(basePath, selectId+settings['SPECS_EXT'])
    selectOutFilePath = \
            os.path.join(basePath, selectId+settings['SELECT_EXT'])

    #print(selectOutFilePath)
    #print(selectId)

    tracesOutFile = TraceFile(tracesOutFilePath)
    specsOutFile = SpecsFile(specsOutFilePath)
    #selectOutFile = SpecsFile(selectOutFilePath)

    tracesOutFile.clear()
    specsOutFile.clear()
    #selectOutFile.clear()

    for mFile in movieFiles:
        trc, spc, _ = get_traces_subset(mFile, settings, selectId)

        if len(trc) > 0:
            tracesOutFile.append(*trc)
            specsOutFile.append(spc)
            #selectOutFile.append(sel)
