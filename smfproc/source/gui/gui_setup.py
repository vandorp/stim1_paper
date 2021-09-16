""" Assemble the GUI """

# System packages
import os
import threading


# Program modules - local relative imports
from . import Tk
from .movies_gui import MoviePanel, MoleculeOverlayPanel
from .traces_gui import TracePanel, StateOverlayPanel
from .core_gui import TkItems, TopPanel, get_monitor_size
from .batch_gui import BatchPanel
from .output_gui import MovieAnalysisPanel



def build_gui():
    main = MainPanel()
    main.set_position([10,30])
    main.title('View')
    
    # === VIEWER PANELS === #

    # Batch panel
    batchPanel = TopPanel(main, 'movies')    
    batchGUI = BatchPanel(batchPanel)
    
    # Movie panel
    moviePanel = TopPanel(main, 'frames')
    
    movieGUI = MoviePanel(moviePanel)
    movieGUI.figurePosition = [830,340]
    movieGUI.figureSize = [6,6]

    moleculeOverlayGUI = MoleculeOverlayPanel(moviePanel, movieGUI)

    # Traces panel
    tracePanel = TopPanel(main, 'traces')

    traceGUI = TracePanel(tracePanel)
    traceGUI.figurePosition = [10,340]
    traceGUI.figureSize = [10,6]

    stateGUI = StateOverlayPanel(tracePanel)
    traceGUI.get_item('listbox').add_update_callback(stateGUI.update)

    # Analysis panel
    analysisPanel = TopPanel(main, 'analysis')
    movieAnalysisGUI = MovieAnalysisPanel(analysisPanel)
    #traceAnalysisGUI = TraceAnalysisPanel(analysisPanel)

    # Output panel
    #outputPanel = TopPanel(main, 'output')
    #outputGUI = OutputPanel(outputPanel)

    # === CALLBACKS & POINTERS === #
    movieList = batchGUI.get_item('movie_list')
    traceList = traceGUI.get_item('listbox')
    
    #overlayGUI.get_subset = traceList.get_subset
    moleculeOverlayGUI.tracePanel = traceGUI

    movieList.add_update_callback(main.reload_file_data)

    traceList.add_update_callback(moleculeOverlayGUI.update)
    traceGUI.add_update_callback(moleculeOverlayGUI.toggle_molecules)

    #batchGUI.traceGUI = traceGUI
    #batchGUI.movieGUI = movieGUI
    #batchGUI.overlayGUI = moleculeOverlayGUI

    #movieGUI.figureHandle.set_size([9,9])
    main.pack()
        
    #traceList.update()

    #main.mainloop()
    return main



class MainPanel(Tk.Tk, TkItems):

    def __init__(self):
        Tk.Tk.__init__(self)
        TkItems.__init__(self)
        
        self.quitButton = Tk.Button(
                master=self, 
                text='Quit', 
                command=self._quit)
        self.bind('q', self._quit)
        
        """
        self.analysisButton = Tk.Button(
                master=self, 
                text='Analysis window', 
                command=self.toggle_analysis_window)
        """

        #self.windows = []
        self.settings = {}
        self.updateRequired = False
        self.threads = []

    def add_thread(self, target, args, descr=''):
        thread = threading.Thread(target=target, args=args)
        self.threads.append([thread, descr])
    
    def execute_threads(self):
        #join_threads = lambda:[thread.join() for thread in self.threads]
        threading.Thread(target=self.join_threads).start()
        self.threads = []

    def join_threads(self):
        self.start_polling()
        for thread in self.threads:
            job, descr = thread
            print('\n=== started %s ===\n'%descr)
            job.start()
            job.join()
            print('\n=== %s finished ===\n'%descr)
        self.stop_polling()

    def start_polling(self):
        pass

    def stop_polling(self):
        pass
    
    """
    def start_polling(self, interval=500):
        self.updateRequired = False
        self.polling = self.after(interval, self.polling_callback, (interval))
    
    def polling_callback(self, args):
        nrThreads = len(threading.enumerate())
        interval = args

        if nrThreads > 1:
            self.polling = self.after(interval, self.polling_callback, (interval))
        #else:
        #    self.after_cancel(self.polling)
        movieList = self.get_item('movies').get_item('batch_sub').get_item('movie_list')
        currMovie = movieList.get()
        if self.updateRequired:
            self.updateRequired = False
            #movieList.callback(currMovie)
            self.reload_file_data(currMovie)
            print 'GUI updated\n'
    """

    def reload_file_data(self, *args):
        batchGUI = self.get_item('movies').get_item('batch_sub')
        traceGUI = self.get_item('traces').get_item('traces')
        moleculeOverlayGUI = self.get_item('frames').get_item('overlay')
        movieGUI = self.get_item('frames').get_item('movie')
        stateOverlayGUI = self.get_item('traces').get_item('state_overlay')
        try:
            fileName = args[0]
            flag = os.path.splitext(fileName)[-1]
        except:
            flag = None

        if flag in self.settings['movie']['EXTENSIONS']:
            #stateOverlayGUI.clear()
            moleculeOverlayGUI.clear()
            traceGUI.reset()
            movieGUI.load_movie(fileName)
            traceGUI.load_movie_traces(fileName)
            moleculeOverlayGUI.load_coordinates(fileName)
            moleculeOverlayGUI.toggle_molecules()
            #stateOverlayGUI.update()

        elif flag==self.settings['core']['TRACE_EXT'] or flag=='.hdf5':
            #movieList = batchGUI.get_item('movie_list')
            #movieList.clear()
            
            movieGUI.get_var('show_movie').set("Skip movie")
            #stateOverlayGUI.clear()
            moleculeOverlayGUI.clear()
            movieGUI.reset()
            movieGUI.close_window()

            traceGUI.reset()
            if flag=='.hdf5':
                traceGUI.load_hdf5(fileName)
            else:
                traceGUI.load_trace_file(fileName)
            #stateOverlayGUI.update()

        else:
            moleculeOverlayGUI.clear()
            movieGUI.reset()
            movieGUI.close_window()
            traceGUI.reset()
                
        
        #print "\nGUI updated"

    def pack(self):
        [p.pack() for p in self.items]
        self.quitButton.pack(side=Tk.BOTTOM, fill=Tk.BOTH)
        #self.analysisButton.pack(side=Tk.BOTTOM, fill=Tk.BOTH)

    def set_position(self, crds):
        #w, h = get_monitor_size()
        self.geometry('+%d+%d'%tuple(crds))
    
        """
        def toggle_analysis_window(self):
        for window in self.windows:
            if window.hidden:
                window.show()
            else:
                window.hide()
        """

    def _quit(self, *args, **kwargs):
        [item._quit() for item in self.items]
        self.quit()
        self.destroy()

