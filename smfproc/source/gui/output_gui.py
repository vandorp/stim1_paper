""" GUI panel for analysis output. """

# System packages
import os, time
import numpy as np
from functools import partial

# Program modules - local relative imports
from . import Tk
from .core_gui import EntryPanel, NumEntry, BoolEntry

# Program modules
from smfproc.source.batch.smfbatch import run_batch_analysis
from smfproc.source.io.filetypes import SpecsFile
from smfproc.source.io.datatypes import SmfSpecs


class MovieAnalysisPanel(EntryPanel):
    """ Controls for movie analysis. """

    def __init__(self, master):
        EntryPanel.__init__(self, master, 'movie')

        cascade = BoolEntry(
                master=self, 
                label="Cascade: ")
        cascade.select()

        distributed = BoolEntry(
                master=self, 
                label="Distributed processing: ")
        distributed.select()

        purge = BoolEntry(
                master=self, 
                label="Purge existing analysis files: ")
        purge.deselect()
        
        manual = BoolEntry(
                master=self, 
                label="Manual selection: ")
        manual.deselect()

        
        registerButton = Tk.Button(
                master=self,
                text='Compile registration files',
                command=partial(self.run_analysis, 'register'))

        laserButton = Tk.Button(
                master=self,
                text='Generate laser profile',
                command=partial(self.run_analysis, 'laser'))

        moviesButton = Tk.Button(
                master=self, 
                text='Analyze movies', 
                command=partial(self.run_analysis, 'movies'))

        tracesButton = Tk.Button(
                master=self, 
                text='Analyze traces', 
                command=partial(self.run_analysis, 'traces'))

        statesButton = Tk.Button(
                master=self, 
                text='Analyze states', 
                command=partial(self.run_analysis, 'states'))
        
        outputButton = Tk.Button(
                master=self, 
                text='Output data & figures', 
                command=partial(self.run_analysis, 'output'))

        
        self.add_item(cascade, 'cascade')
        self.add_item(distributed, 'distributed')
        self.add_item(purge, 'purge')
        self.add_item(manual, 'manual')

        self.add_item(registerButton, 'register')
        self.add_item(laserButton, 'laser')
        self.add_item(moviesButton, 'movies')
        self.add_item(tracesButton, 'traces')
        self.add_item(statesButton, 'states')
        self.add_item(outputButton, 'output')

        self.main = master.master

    def get_movie_list(self):
        movieGUI = self.main.get_item('frames').get_item('movie')
        batchGUI = self.main.get_item('movies').get_item('batch_sub')
        movieList = batchGUI.get_item('movie_list')
        return movieList.content, movieGUI

    def run_analysis(self, descr):
        moviePaths, movieGUI = self.get_movie_list()

        cascade = self.get_item('cascade').get()
        purge = self.get_item('purge').get()
        distributed = self.get_item('distributed').get()
        settings = self.main.settings
        manual = self.get_item('manual').get()

        if distributed:
            nr_cores = 0
        else:
            nr_cores = 1

        func = run_batch_analysis
        args = (moviePaths, descr, purge, False, cascade, nr_cores, manual, settings)

        if descr=='output':
            # Run in main thread
            func(*args)
        else:
            # Run in separate thread
            self.main.add_thread(target=func, args=args, descr='batch analysis')
            self.main.execute_threads()
            self.main.start_polling()

