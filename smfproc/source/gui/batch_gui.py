""" BATCH GUI """

# System packages
import os, glob
import numpy as np
import threading

# Program modules - local relative imports
from . import Tk, tkFileDialog
from .core_gui import BoolEntry, NumEntry, SubPanel, ListBox

# Program modules
from smfproc.source.io.filetypes import TraceFile, SpecsFile
from smfproc.source.io.misc import load_file_paths
from smfproc.source.analysis.movies.smfmovie import analyze as analyze_movie
from smfproc.source.batch.batch_misc import export_traces_subset, get_shared_base_path


class BatchPanel(SubPanel):
    """ GUI panel for loading movies and traces. """

    def __init__(self, master):
        SubPanel.__init__(self, master, descr='batch_sub')

        
        movieList = ListBox(master=self, width=30, height=8)
        #traceList = ListBox(master=self)
        
        """
        traceDataFileVar = Tk.StringVar()
        traceDataFile = Tk.Entry(
                master=self,
                textvariable=traceDataFile)
        """

        #movieList.add_update_callback(self.update)

        #self.movieExtensions = []
        #self.sources = {}
        self.main = master.master
        #self.traceGUI = None
        #self.movieGUI = None
        #self.overlayGUI = None
        
        loadBatchButton = Tk.Button(
                master=self, 
                text='Load movie batch', 
                command=self.load_batch)

        loadFileButton = Tk.Button(
                master=self, 
                text='Load movie/trace file', 
                command=self.load_file)

        clearListButton = Tk.Button(
                master=self,
                text='Clear list',
                command=self.clear)

        #exportTracesButton = Tk.Button(
        #        master=self, 
        #        text='Export current trace', 
        #        command=self.export_traces)

        self.add_item(movieList, 'movie_list')
        #self.add_item(traceList, 'trace_list')
        #self.add_item(traceDataFile, 'data_file', traceDataFileVar)
        self.add_item(loadBatchButton, 'load_batch')
        self.add_item(loadFileButton, 'load_file')
        self.add_item(clearListButton, 'clear_list')
        #self.add_item(exportTracesButton, 'export_traces')
       
        #self.update(moviePaths)

    def export_traces(self):
        self.traceGUI = self.main.get_item('traces').get_item('traces')
        """
        self.traceGUI.export_selections()

        selectId = self.traceGUI.get_var('select').get()
        #selectId = 'selected'
        
        movieFiles = self.get_item('movie_list').content
        """
        fileOutName = self.get_item('movie_list').listBox.get(0)
        self.traceGUI.export_current_trace(fileOutName)
        #print(os.curdir)

        #sharedPath = get_shared_base_path(movieFiles)
        #basePath = os.path.sep.join(sharedPath[:-1])

        #settings = self.main.settings['core']
        
        #export_traces_subset(selectId, movieFiles, settings)

        """
        tracesOutFilePath = \
                os.path.join(basePath, selectId+settings['TRACE_EXT'])
        specsOutFilePath = \
                os.path.join(basePath, selectId+settings['SPECS_EXT'])
        selectOutFilePath = \
                os.path.join(basePath, selectId+settings['SELECT_EXT'])

        tracesOutFile = TraceFile(tracesOutFilePath)
        specsOutFile = SpecsFile(specsOutFilePath)
        selectOutFile = SpecsFile(selectOutFilePath)
        
        tracesOutFile.clear()
        specsOutFile.clear()
        selectOutFile.clear()

        for mFile in movieFiles:
            trc, spc, _ = get_traces_subset(mFile, settings, selectId)
            
            tracesOutFile.append(*trc)
            specsOutFile.append(spc)
            #selectOutFile.append(sel)
        """

    def update(self, content=list()):
        movieList = self.get_item('movie_list')

        if len(content) > 1:
            sharedPath = os.path.sep.join(get_shared_base_path(content)[:-1])
            #descriptions = [fPath.split(sharedPath)[1].split(os.path.sep)[1] 
            #        for fPath in content]

            descriptions = [fPath.split(sharedPath)[1] for fPath in content]
            #descriptions = [os.path.splitext(os.path.basename(fPath))[0] 
            #        for fPath in content]

            allSameFiles = sum([len(descr)>0 for descr in descriptions])==0

            if allSameFiles:
                fileName = os.path.basename(content[0])
                descriptions = [fileName]*len(content)
        else:
            descriptions = [fPath.split(os.path.sep)[-1] for fPath in content]

        movieList.set(content=content, display=descriptions)
        movieList.show()
        # NB movie doesn't refresh here, but don't understand why not.
        # This is a little hack for now.
        if len(content) > 0:
            #    movieList.callback(movieList.content[0])
            #movieList.listBox.select_set(len(content)-1) #This only sets focus on the first item.

            #movieList.listBox.event_generate("<<ListboxSelect>>")
            movieList.listBox.selection_clear(0)
            movieList.listBox.select_set(0) #This only sets focus on the first item.
            self.main.reload_file_data(movieList.content[0])

        # Set focus on the list of traces
        traceList = self.main.get_item('traces').get_item('traces').get_item('listbox')
        traceList.take_focus()


    def load_batch(self):
        dirPath = tkFileDialog.askdirectory()
        if len(dirPath)==0:
            return

        movieExtensions = self.main.settings['movie']['EXTENSIONS']
        moviePaths = load_file_paths(dirPath, None, movieExtensions)

        movieList = self.get_item('movie_list')
        currentMoviePaths = movieList.content

        newMovieList = moviePaths + currentMoviePaths
        self.update(newMovieList)

    def load_file(self):
        #self.traceGUI = self.main.get_item('traces').get_item('traces')
        #self.moleculeOverlayGUI = self.main.get_item('frames').get_item('overlay')
        #self.movieGUI = self.main.get_item('frames').get_item('movie')
        #self.stateOverlayGUI = self.main.get_item('traces').get_item('state_overlay')

        fileName = tkFileDialog.askopenfilename()

        if len(fileName) == 0:
            return

        #self.main.reload_file_data(fileName)

        movieList = self.get_item('movie_list')
        
        flag = os.path.splitext(fileName)[-1]
        #if flag in self.main.settings['movie']['EXTENSIONS']:
        currentMoviePaths = movieList.content
        newMovieList = [fileName] + currentMoviePaths
        self.update(newMovieList)

        #elif flag=='.trc':
        #    self.update([fileName])

            
        #else:
        #    print 'Error: unknown file type.'

    def clear(self):
        self.update([])
        self.main.reload_file_data()
        

