""" Builds a GUI panel for movie viewing. """

# System packages
import os, shutil
import numpy as np
import time
import threading
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Ellipse
from matplotlib import gridspec, rcParams

# Program modules
from smfproc.source.batch.batch_analysis import load_movie_data
from smfproc.source.analysis.movies.analysis import get_max_fluorescence
from smfproc.source.io.filetypes import SpecsFile, TraceFile
from smfproc.source.analysis.movies.smfmovie import analyze

# Program modules - local relative imports
from . import Tk
from .core_gui import NumEntry, BoolEntry, TopPanel, SubPanel, EntryPanel, \
        FigurePanel, ScrollBarGUI, OptionMenu
from .core_figures import FigureBase




rcParams.update({'font.size': 12})


class MoviePanel(FigurePanel):
    """ Builds the GUI panel for viewing smFRET movies. """

    def __init__(self, master):
        FigurePanel.__init__(self, master, descr='movie')
        
        # === Create panel items
        scrollBar = ScrollBarGUI(master=self)
        scrollBar.add_update_callback(self.update)
 
        metaDataDisplay = Tk.Text(
                master=self,
                width=20,
                height=10)
        
        """
        toggleMovieVar = Tk.StringVar()
        toggleMovie = Tk.Checkbutton(
                master=self, 
                text="Show movies",
                variable=toggleMovieVar)
        toggleMovie.select()
        """
        
        #loadLabel = Tk.Label(master=self, 
        #        text='Load: ',
        #        anchor=Tk.W)
        
        selectMovieVar = Tk.StringVar()
        selectMovie = OptionMenu(
                self,
                selectMovieVar,
                "Skip movie", "First frames only", "Load full movie", "Maximum signal")
        #selectMenuVar.set("")
        #selectMenuVar.trace('w', self.callback)

        self.add_item(scrollBar, 'scrollbar')
        self.add_item(metaDataDisplay, 'display')
        #self.add_item(toggleMovie, 'show_movie', toggleMovieVar)
        #self.add_item(loadLabel, 'load_label')
        self.add_item(selectMovie, 'show_movie', selectMovieVar)

        # Open a window here to allow overlays
        #self.open_window(ImageFigure)
        
        #self.figureHandle.show()
        #self.donAccOverlay = OverlayFigure(self.figureHandle)
        
        #self.moleculeOverlays = []

        # === Keyboard bindings
        #master.bind('<Left>', scrollBar.key_prev)
        #master.bind('<Right>', scrollBar.key_next)

        #self.sources = {}
        #self.active = True

        self.main = self.master.master
        self.currentMoviePath = None

        """
        def toggle_movies(self):
        #movieList = self.get_item('movie_list')
        #if len(movieList.content) > 0:
        #    movieList.callback(movieList.content[0])
        self.active = self.get_var("show_movie").get()=='0'
        """
    def is_active(self):
        return not self.get_var("show_movie").get()=='Skip movie'

    def reset(self):
        if self.has_data():
            self.data = []

        if self.has_window():
            self.figureHandle.blank()
            self.figureHandle.set_title('')
            self.figureHandle.show()
        
        scroll = self.get_item('scrollbar')
        scroll.config(0)
        self.reset_meta_display()

    def load_movie(self, filePath):
        
        if not self.is_active():
            self.close_window()
            return

        #if self.currentMoviePath == filePath:
        #    return

        # Reset the GUI
        self.reset()
        
        #if self.has_data():
        self.reset_meta_display('Movie loading...')
        #self.figureHandle.set_title('Movie loading...')

        settings = self.main.settings['core']
        
        metaFile = settings['FILE_BASE']+settings['META_EXT']
        metaPath = os.path.join(os.path.splitext(filePath)[0], metaFile)
        movieMeta = TraceFile(metaPath).read()

        # Load movie data
        movieMode = self.get_var("show_movie").get() 
        if movieMode == "First frames only":
            #nrFrames = 5
            movieData, movieHeader = load_movie_data(filePath, 5)
        elif movieMode == "Load full movie":
            #nrFrames = -1
            movieData, movieHeader = load_movie_data(filePath, -1)
        elif movieMode == "Maximum signal":
            movs, movieHeader = load_movie_data(filePath, -1)

            don = movs[0]
            acc = movs[1]

            donMaxIm = get_max_fluorescence(don.frames, 5)
            accMaxIm = get_max_fluorescence(acc.frames, 5)

            don.frames = np.reshape(donMaxIm, [donMaxIm.shape[0], donMaxIm.shape[1], 1])
            acc.frames = np.reshape(accMaxIm, [accMaxIm.shape[0], accMaxIm.shape[1], 1])

            movieData = [don, acc]
        else:
            return


        #guiMode = self.main.get_item('analysis').get_item('movie').get_item('register').cget('text')
        #if guiMode == 'Registration mode':
        #    self.get_var("show_movie").set("Load full movie")
        #    nrFrames = -1
        
        #movieMeta = load_movie_meta(filePath)

        if not self.has_window():
            self.open_window(ImageFigure)

        nrFrames = movieData[0].nr_frames()
        metaData = self.compile_meta_data(nrFrames, movieHeader, movieMeta)

        self.data = [movieData, metaData]
        self.currentMoviePath = filePath

        if self.has_data():
            scroll = self.get_item('scrollbar')
            scroll.config(nrFrames)
            #display.delete(1.0, Tk.END)

            ind = scroll.get()
            self.update(val=ind, redraw=True)
        else:
            self.reset()
            #self.close_window()
            #scroll.config(0)

    def update(self, val, redraw=False):
        if not self.is_active():
            return

        if ((not self.has_window()) | (not self.has_data())):
            return

        movie = self.data[0]
        
        nrFrames = movie[0].nr_frames()

        frames = [m.get_frame(val) for m in movie]

        if redraw:
            self.figureHandle.cla()
            self.figureHandle.draw(frames)
        else:
            self.figureHandle.update(frames)
        
        self.figureHandle.set_title('frame %d of %d' % (val+1, nrFrames))
        self.figureHandle.show()
     
        self.update_meta_display(val)

    def compile_meta_data(self, nrFrames, perMovie={}, perFrame=[]):
        meta = {}

        for key in perMovie:
            meta[key] = [str(perMovie[key])]*nrFrames
        
        #if len(perFrame) > 0:
        laserTs = perFrame[0]
        laserDescr = perFrame[1]
        for n in range(len(laserTs)):
            meta[laserDescr[n]] = [str(p) for p in laserTs[n].data[0]]

        return meta
 
    def update_meta_display(self, frameIndex=0):
        display = self.get_item('display')
        display.config(state=Tk.NORMAL)
        display.delete("1.0", Tk.END)
        meta = self.data[1]
        for key in sorted(meta.keys()):
            display.insert(Tk.END, key+': '+meta[key][frameIndex]+'\n')
        display.config(state=Tk.DISABLED)

    def reset_meta_display(self, content=''):
        display = self.get_item('display')
        display.config(state=Tk.NORMAL)
        display.delete("1.0", Tk.END)
        display.insert(Tk.END, content+'\n')
        display.config(state=Tk.DISABLED)
        display.update()




class MoleculeOverlayPanel(EntryPanel):
    """ Makes a movie frame overlay for indicating single molecules. """

    def __init__(self, master, moviePanel):
        SubPanel.__init__(self, master, descr='overlay')
        
        """
        toggleMoleculesVar = Tk.StringVar()
        toggleMolecules = Tk.Checkbutton(
                master=self, 
                text="Show molecules",
                variable=toggleMoleculesVar,
                command=self.toggle_molecules)
        toggleMolecules.select()
        """
        toggleMolecules = BoolEntry(
                master=self,
                label="Show molecules",
                command=self.toggle_molecules)
        toggleMolecules.select()

        #toggleMolecules.config(state=Tk.DISABLED)
        #self.molecule_toggle_callbacks = []
        
        self.add_item(toggleMolecules, 'molecules')
        
        self.moviePanel = moviePanel
        #self.figureHandle = moviePanel.figureHandle
        #self.donorOverlay = MoleculeOverlay(self.figureHandle.axes[0])
        #self.acceptorOverlay = MoleculeOverlay(self.figureHandle.axes[1])
        #self.donorOverlay = MoleculeOverlay()
        #self.acceptorOverlay = MoleculeOverlay()

        #self.sourceFile = None
        #self.sources = {}
        self.tracePanel = None # Pointer to tracePanel
        
        #def get_subset(self):
        #return np.array([])

        self.main = self.master.master

        """
        def show_subset(self, *args):
        print self.tracePanel.get_item('listbox').subset
        #print self.get_subset()
        """

    def toggle_molecules(self, *args, **kwargs):
        if not self.moviePanel.is_active():
            return

        if not self.moviePanel.has_window():
            return

        toggleState = self.get_item('molecules').get()
        #for ovl in self.molecule_overlays:
        #    ovl.toggle_show_molecules(state)

        subsetState = self.tracePanel.get_item('listbox').subset == 1
        
        fig = self.moviePanel.figureHandle

        if len(subsetState) == len(fig.donorOverlay.overlay):
            state = subsetState*toggleState
        else:
            state = toggleState

        fig.donorOverlay.toggle(state)
        fig.acceptorOverlay.toggle(state)

        #self.update()

        fig.show()

    def load_coordinates(self, moviePath):
        if not self.moviePanel.is_active():
            return
        
        if not self.moviePanel.has_window():
            return

        settings = self.main.settings['core']
        
        specsFile = settings['FILE_BASE']+settings['SPECS_EXT']
        specsPath = os.path.join(os.path.splitext(moviePath)[0], 
                specsFile)
        specs = SpecsFile(specsPath).read()

        if len(specs.descr) == 0:
            return
        
        fig = self.moviePanel.figureHandle
        #specs = load_trace_specs(filePath)
        fig.donorOverlay.clear()
        fig.donorOverlay.draw(specs.get('don_row'), specs.get('don_col'), 
                specs.get('radius')+1)
        fig.acceptorOverlay.clear()
        fig.acceptorOverlay.draw(specs.get('acc_row'), specs.get('acc_col'), 
                specs.get('radius')+1)
        
        #self.show_subset()
        self.update()
        fig.show()

    def update(self, *args, **kwargs):
        if not self.moviePanel.is_active():
            return

        if not self.moviePanel.has_window():
            return
        
        fig = self.moviePanel.figureHandle

        index = self.tracePanel.get_item('listbox').get()
        fig.donorOverlay.update(-1, edgecolor='y', linewidth=1)
        fig.donorOverlay.update(index, edgecolor='g', linewidth=2)
        
        fig.acceptorOverlay.update(-1, edgecolor='y', linewidth=1)
        fig.acceptorOverlay.update(index, edgecolor='r', linewidth=2)

        fig.show()

    def clear(self, *args, **kwargs):
        if self.moviePanel.has_window():
            fig = self.moviePanel.figureHandle
            fig.donorOverlay.clear()
            fig.acceptorOverlay.clear()


class ImageFigure(FigureBase):
    """ Class for showing figure with movie frame image. """

    def __init__(self, figsize=(6,6)):
        figsize = [6,6]
        FigureBase.__init__(self, figsize=figsize)
        #self.figure, (self.axDon, self.axAcc) = \
        #        plt.subplots(1, 2, sharex=True, sharey=True, figsize=figsize)
        
        #self.figure = plt.figure(figsize=figsize)
        pos = gridspec.GridSpec(1,2)
        
        self.axes = []
        self.axes.append(self.figure.add_subplot(pos[0]))
        self.axes.append(self.figure.add_subplot(pos[1], 
            sharex=self.axes[0], sharey=self.axes[0]))

        self.donorOverlay = MoleculeOverlay(self.axes[0])
        self.acceptorOverlay = MoleculeOverlay(self.axes[1])


    def draw(self, ims): 
        # Show images
        self.imDon = self._show_image(self.axes[0], ims[0])
        self.imAcc = self._show_image(self.axes[1], ims[1])

        # Report image values
        get_pixel_value = lambda x, y, im: im.get_array()[int(x), int(y)]

        self.axes[0].format_coord = lambda x, y: \
                "row: {0:0.0f}  ".format(y) + \
                "col: {0:0.0f}  ".format(x) + \
                "value: {0:0.0f}".format(get_pixel_value(y, x, self.imDon))
        self.axes[1].format_coord = lambda x, y: \
                "row: {0:0.0f}  ".format(y) + \
                "col: {0:0.0f}  ".format(x) + \
                "value: {0:0.0f}".format(get_pixel_value(y, x, self.imAcc))
        
        # Figure layout
        self.set_layout()
        self.figure.subplots_adjust(0.01,0.01,0.98,0.98,0,0)
        
    @staticmethod
    def _show_image(ax, im):
        ax.set_axis_off()
        plt.setp(ax, aspect=1.0, adjustable='box-forced')
        img = im.get()
        return ax.imshow(img, cmap='jet')
        #return ax.imshow(img, cmap='gray')

    def set_layout(self):
        data = self.imDon.get_array()
        self.axes[0].axis([0, data.shape[1], data.shape[0], 0])
        data = self.imAcc.get_array()
        self.axes[1].axis([0, data.shape[1], data.shape[0], 0])

    def update(self, ims):
        self.imDon.set_data(ims[0].get())
        self.imAcc.set_data(ims[1].get())

    def cla(self):
        [ax.cla() for ax in self.axes]

    def blank(self):
        don = self.imDon.get_array()
        acc = self.imAcc.get_array()
        self.imDon.set_data(don*0)
        self.imAcc.set_data(acc*0)


class MoleculeOverlay():
    """ Overlay for movie frame to indicate single molecules. """

    def __init__(self, ax=None):
        self.ax = ax
        self.overlay = []

    def draw(self, rows, cols, radii):
        #for r, c in zip(crds[0], crds[1]):
        #    mol = self._draw_molecule(c, r, radius)
        #    self.overlay.append(mol)

        for n in range(len(rows)):
            mol = self._draw_molecule(cols[n], rows[n], radii[n])
            self.overlay.append(mol)

        #self.leftHighlight = ImageFigure._draw_molecule(self.axes[0], 
        #        rcLeft[1][0], rcLeft[0][0], radius, (0,1,0), linewidth=2)
       
    def _draw_molecule(self, col, row, radius, color='yellow', *args, **kwargs):
        return self.ax.add_artist(Circle(xy=(col, row), 
            radius=radius, 
            edgecolor=color,
            facecolor='none', 
            visible=True,
            *args, **kwargs))

    def clear(self):
        try:
            [mol.remove() for mol in self.overlay]
        except:
            pass
        self.overlay = []

    def update(self, index, *args, **kwargs):
        if not index==None:
            if index==-1:
                [m.set(*args, **kwargs) for m in self.overlay]
            else:
                self.overlay[index].set(*args, **kwargs)

    def toggle(self, visible):
        if type(visible) == bool:
            visible = visible*np.ones(len(self.overlay))
        [self.overlay[m].set_visible(visible[m]) for m in range(len(self.overlay))]

