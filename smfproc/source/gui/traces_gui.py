""" GUI for viewing traces. """


# System packages
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
from pylab import get_current_fig_manager

# Program modules - local relative imports
from . import Tk
from .core_gui import SubPanel, EntryPanel, NumEntry, FigurePanel, \
        SuperListBox, OptionMenu, BoolEntry
from .core_figures import FigureBase
from .drag_resize_rectangle import DraggableResizeableRectangle as drag

# Program modules
from smfproc.source.io.filetypes import SpecsFile, TraceFile
from smfproc.source.io.datatypes import State, SmfSpecs
from smfproc.source.analysis.traces.analysis import correct_donor, correct_fret
from smfproc.source.analysis.traces.smftraces import collect_trace_data, analyze as analyze_traces
from smfproc.source.analysis.states.analysis import get_state_model


rcParams.update({'font.size': 12})



class TracePanel(FigurePanel):
    """ GUI panel for viewing traces. """

    def __init__(self, master):
        FigurePanel.__init__(self, master, descr='traces')
        
        # === Create panel items
        listBox = SuperListBox(master=self)
        listBox.add_update_callback(self.update)
        listBox.superSelectedCallback = self.export_selections

        selectMenuVar = Tk.StringVar()
        selectMenu = OptionMenu(
                self,
                selectMenuVar,
                "")
        #selectMenuVar.set("")
        selectMenuVar.trace('w', self.callback)

        registerButton = Tk.Button(
                master=self,
                text='Export registration file',
                command=self.export_registration)
        
        exportTracesButton = Tk.Button(
                master=self, 
                text='Export current trace', 
                command=self.export_current_trace)

        self.add_update_callback(self.update_selection)

        self.add_item(listBox, 'listbox')
        self.add_item(selectMenu, 'select', selectMenuVar)
        #self.add_item(exportButton, 'export')
        self.add_item(registerButton, 'register')
        self.add_item(exportTracesButton, 'export_traces')

        self.filePath = ''
        #self.selections = SmfSpecs()
        #self.sourceFile = ''
        #self.sources = {}
        
        self.master = master
        self.main = master.master

    def _quit(self):
        FigurePanel._quit(self)
        self.export_selections()
    
    """
    def export_traces(self):
        listBox = self.get_item('listbox')
        subset = listBox.subset

        traces = [d.select(subset) for d in self.data[0]]
        
        selected = self.get_var('select').get()
        outFileName = selected + \
                self.main.settings['core']['TRACE_EXT']

        tracePath = os.path.join(os.path.split(self.filePath)[0], 
                outFileName)

        outFile = TraceFile(tracePath)
        #outFile.write(traces, {'blocks':self.data[1]})
        outFile.write(traces, self.data[1])
    """

    def export_current_trace(self):
        #print self.filePath
        moviePanel = self.main.get_item('movies').get_item('batch_sub')
        fileOutBase = moviePanel.get_item('movie_list').listBox.get(0)
        #fileOutBase = os.path.dirname(self.filePath)

        listBox = self.get_item('listbox')
        val = listBox.get()

        traces = [ts.get(val) for ts in self.data[0]]
        descr = self.data[1]
        
        fileOutName = os.path.splitext(fileOutBase)[0] + '_molecule_%d'%(val+1) + '.trc'
        fileOutName = fileOutName.replace(os.path.sep, '_')
        tracesOutFile = TraceFile(fileOutName)
        tracesOutFile.write(traces, descr)
    
    def export_selections(self):
        listBox = self.get_item('listbox')
        
        selections = self.data[-1]
        #selections.set(np.ones(listBox.size()), 'all')
        selections.set(listBox.get_super_selected(), 'manual')
        #selections.update_selected()
        
        finalSelected = selections.get(['manual', 'auto'], np.logical_xor)
        selections.set(finalSelected, 'selected')

        selectPath = os.path.splitext(self.filePath)[0] + \
                self.main.settings['core']['SELECT_EXT']

        SpecsFile(selectPath).write(selections)
    
    def export_registration(self):

        #specs = self.main.get_item('traces').get_item('traces').data[-2]
        #selections = self.main.get_item('traces').get_item('traces').data[-1]
        specs = self.data[-2]
        selections = self.data[-1]

        #cat = self.main.get_item('traces').get_item('traces').get_var('select').get()
        cat = self.get_var('select').get()

        manual = selections.get('manual')
        regist = selections.get(cat)
        store = np.logical_and(manual==False, regist)

        don_row = specs.get('don_row')[store]
        don_col = specs.get('don_col')[store]
        acc_row = specs.get('acc_row')[store]
        acc_col = specs.get('acc_col')[store]

        #filePath = self.main.get_item('traces').get_item('traces').filePath
        filePath = self.filePath
        outPath = os.path.join(os.path.dirname(filePath), 'registration.specs')

        #print outPath
        spc = SmfSpecs()
        spc.set(don_row, 'don_row')
        spc.set(don_col, 'don_col')
        spc.set(acc_row, 'acc_row')
        spc.set(acc_col, 'acc_col')

        SpecsFile(outPath).write(spc)

        print('Registration written as %s\n'%outPath)
    
    def reset(self, *args, **kwargs):
        if self.has_window():
            stateOverlayGUI = self.master.get_item('state_overlay')

            stateOverlayGUI.clear()

            self.figureHandle.cla()
            self.figureHandle.set_layout()
            self.figureHandle.set_title('')
            self.figureHandle.show()

        if self.has_data():
            #self.export_selections()
            self.data = []
            self.get_item('listbox').clear()
            self.get_item('select').clear()
            self.get_var('select').set("")    
    
    def load_hdf5(self, filePath):
        #need python3
        pass

    def load_trace_file(self, filePath):
        
        self.filePath = filePath
                
        settings = self.main.settings['core']

        traces, descr, specs, selections = \
                collect_trace_data(filePath, settings)

        if len(traces) == 0:
            return

        gamma = specs.get('gamma')
        
        donCorr = traces[0].copy()
        correct_donor(donCorr.data, gamma)

        fretCorr = traces[4].copy()
        correct_fret(fretCorr.data, gamma)

        traces.extend([donCorr, fretCorr])
        descr.extend(['donor_corrected', 'fret_corrected'])
        #traces = [donCorr, fretCorr] + traces
        #descr = ['donor_corrected', 'fret_corrected'] + descr

        data = [traces, descr, specs, selections]

        self._fill_gui(data)
        
        # Fill data last to prevent running self.update_selection too early
        self.data = data
        self.update_selection()

    def load_movie_traces(self, moviePath):
        baseDir = os.path.splitext(moviePath)[0]

        settings = self.main.settings['core']
        
        tracePath = os.path.join(baseDir, 
                settings['FILE_BASE']+settings['TRACE_EXT'])
        
        self.load_trace_file(tracePath)

    def _fill_gui(self, data):
        if len(data) == 0:
            return

        traces = data[0]
        specs = data[2]
        selections = data[3]

        if len(traces) == 0:
            return
          
        if not self.has_window():
            self.open_window(TraceFigure)
        
        nrTraces = traces[0].get_nr_traces()
        # --- Set selections
        #selections.prefill(traces[0].get_nr_traces())
        descr = selections.descr
        if not ('manual' in descr):
            selections.set(np.zeros(nrTraces), 'manual')
        if not ('auto' in descr):
            selections.set(np.ones(nrTraces), 'auto')


        # --- Fill selections menu
        menu = self.get_item('select')
        var = self.get_var('select')
         
        menu.clear()
        [menu.add_option(item, var) for item in descr]
        if np.sum(selections.get('selected'))==0:
            var.set('all')
        else:
            var.set('selected')
        
        # --- Fill listbox
        #indices = np.arange(specs.data.shape[0])
        indices = np.arange(specs.get_nr_rows())
        listbox = self.get_item('listbox')
        listBoxEntries = ['Molecule %d'%(index+1) for index in indices]

        listbox.set(content=indices, display=listBoxEntries)
        listbox.set_super_selected(selections.get('manual'))
        listbox.mode = selections.get('auto')
        listbox.show()

        
    def update_selection(self, *args):
        if not self.has_data():
            return

        # Get associated binary array
        selections = self.data[-1]

        # Make subset of indices
        listbox = self.get_item('listbox')
 
        # Update the stored manual selections
        selections.set(listbox.get_super_selected(), 'manual')

        # Get current selection from menu
        selectionDescr = self.get_var('select').get()
        selectionCrit = selections.get(selectionDescr)

        listbox.show(subset=selectionCrit)

        # Update trace plot
        index = listbox.get()
        self.update(index, replot=True)

    def update(self, val, replot=False):
        if not self.has_data():
            return

        nrTraces = self.data[0][0].get_nr_traces()
        
        if val==None:
            self.figureHandle.cla()
            self.figureHandle.set_layout()
            self.figureHandle.set_title('')
        else:
            self.update_plot(val, replot)    
            self.figureHandle.set_title('trace %d of %d' %
                    (val+1, nrTraces))

        self.figureHandle.show()

    def update_plot(self, val, replot):
        traces = [ts.get(val) for ts in self.data[0]]
        descr = self.data[1]
        
        #stateOverlayGUI = self.main.get_item('traces').get_item('state_overlay')
        stateOverlayGUI = self.master.get_item('state_overlay')

        #stateOverlayGUI.clear

        if replot:
            self.figureHandle.cla()
            self.figureHandle.plot(traces, descr)
            self.figureHandle.set_layout()
        else:            
            self.figureHandle.update(traces, descr)
            self.figureHandle.axes[0].scale_yaxis()
        
        stateOverlayGUI.update()
        stateOverlayGUI.toggle_model()

        #print [overlay.rect.get_x() for overlay in self.figureHandle.accBleachOverlay.overlay]


class TraceSubplot():
    """ Trace figure axis. """

    def __init__(self, fig, *args, **kwargs):
        self.ax = fig.add_subplot(*args, **kwargs)

        self.ax.format_coord = lambda x, y: \
                "x: {0:.1f}  ".format(x) + \
                "y: {0:.1f}".format(y)

        self.cla()

    def cla(self):
        self.ax.cla()
        self.plots = {}
        self.markers = {}

    def add_plot(self, handle, descr):
        self.plots[descr] = handle

    def update_plot(self, data, descr):
        handle = self.plots[descr]
        nrTraces = len(handle)
        if nrTraces==1:
            handle[0].set_ydata(data)
        else:
            for n in range(nrTraces):
                handle[n].set_ydata(data[:,n])

    def scale_yaxis(self):
        if len(self.plots.values()) == 0:
            ymax = 1.0
            ymin = 0
        else:
            yData = np.array([p[0].get_ydata() for p in self.plots.values()])
            ymax = np.max(yData)
            ymin = np.min(yData)

            ymax *= 1.2
            if ymin > -0.1*ymax:
                ymin = -0.1*ymax
            else:
                ymin *= 1.2
        
        self.ax.set_ylim((ymin, ymax))

    def scale_xaxis(self):
        if len(self.plots.values()) > 0:
            xmax = self.ax.axes.dataLim.xmax
            xmin = self.ax.axes.dataLim.xmin
            self.ax.set_xlim((xmin, xmax))
        
    def set_layout(self):
        self.scale_xaxis()
        self.scale_yaxis()

        self.ax.plot(self.ax.get_xlim(), [0, 0],
                color='k', alpha=0.25, linestyle='--', zorder=0)
        

class TraceFigure(FigureBase):
    """ Trace figure window. """

    specs = {
            'donor': (0, 'g', 1),
            'donor_corrected': (0, [0.7, 1, 0.7], 0),
            'acceptor': (0, 'r', 2),
            'fret': (1, 'k', 1),
            'fret_corrected':(1, [0.7, 0.7, 0.7], 0),
            'dy': (2, 'orange', 0),
            'dx': (2, 'grey', 1),
            'laser532': (3, 'g', 0),
            'laser638': (3, 'r', 1)}

    def __init__(self, figsize=(9,6)):
        FigureBase.__init__(self, figsize=figsize)

        self.axes = []

        pos = gridspec.GridSpec(4, 1, height_ratios=[2,2,1,1])
        self.axes.append(TraceSubplot(self.figure, pos[0]))
        self.axes.append(TraceSubplot(self.figure, pos[1], sharex=self.axes[0].ax))
        self.axes.append(TraceSubplot(self.figure, pos[2], sharex=self.axes[0].ax))
        self.axes.append(TraceSubplot(self.figure, pos[3], sharex=self.axes[0].ax))
 
        self.donModelOverlay = ModelOverlay(self.axes[0].ax)
        self.accModelOverlay = ModelOverlay(self.axes[0].ax)
        self.fretModelOverlay = ModelOverlay(self.axes[1].ax)
        self.clusModelOverlay = ModelOverlay(self.axes[1].ax)

        self.donBleachOverlay = BleachZoneOverlay(self.axes[0].ax)
        self.accBleachOverlay = BleachZoneOverlay(self.axes[0].ax)

    def set_layout(self):
        self.axes[0].ax.set_ylabel('Intensity (a.u.)')
        self.axes[1].ax.set_ylabel('FRET ratio')
        self.axes[2].ax.set_ylabel('$\Delta$ (pixels)')
        self.axes[3].ax.set_ylabel('P (mW)')
        self.axes[3].ax.set_xlabel('Time (s)')
        
        #[ax.ax.legend(ax.plots.keys(), fontsize=8) for ax in self.axes]
        [ax.ax.legend(fontsize=8) for ax in self.axes]
        [ax.set_layout() for ax in self.axes]

        self.axes[1].ax.set_ylim([-0.25, 1.25])
        
        self.figure.tight_layout()

    def plot(self, traces, descr):
        for n in range(len(descr)):
            #if TraceFigure.specs.has_key(descr[n]):
            if descr[n] in TraceFigure.specs.keys():
                index = TraceFigure.specs[descr[n]][0]
                color = TraceFigure.specs[descr[n]][1]
                zorder = TraceFigure.specs[descr[n]][2]
                self.axes[index].add_plot(
                        traces[n].plot(self.axes[index].ax, 
                        color=color, zorder=zorder, label=descr[n]), 
                        descr[n])

    def update(self, traces, descr):
        for n in range(len(descr)):
            #if TraceFigure.specs.has_key(descr[n]):
            if descr[n] in TraceFigure.specs.keys():
                index = TraceFigure.specs[descr[n]][0]
                #color = TraceFigure.specs[descr[n]][1]
                data = traces[n].data
                data[np.isnan(data)] = 0
                self.axes[index].update_plot(traces[n].data, descr[n])
        
    def cla(self):
        [ax.cla() for ax in self.axes]


class StateOverlayPanel(EntryPanel):
    """ Toggle trace axis overlays. """

    def __init__(self, master):
        EntryPanel.__init__(self, master, descr='state_overlay')

        toggleModel = BoolEntry(
                master=self, 
                label="Show model: ",
                command=self.toggle_model)
        toggleModel.select()
        
        self.add_item(toggleModel, 'model')
        
        self.main = self.master.master
        self.tracePanel = master.get_item('traces')

    def toggle_model(self, *args, **kwargs):
        if not self.tracePanel.has_window():
            return

        toggleState = self.get_item('model').get()

        fig = self.tracePanel.figureHandle

        fig.donModelOverlay.toggle(toggleState)
        fig.accModelOverlay.toggle(toggleState)
        fig.fretModelOverlay.toggle(toggleState)
        fig.clusModelOverlay.toggle(toggleState)

        fig.show()

    def update(self, *args):
        if not self.tracePanel.has_window():
            return

        index = self.tracePanel.get_item('listbox').get()

        traces = self.tracePanel.data[0]
        descr = self.tracePanel.data[1]

        if not ('states' in descr):
            return
        
        stt = State(traces[descr.index('states')])
        
        don = traces[descr.index('donor')]
        acc = traces[descr.index('acceptor')]
        fret = traces[descr.index('fret')]
        
        fig = self.tracePanel.figureHandle
        self.clear()

        # Model overlay
        donSigModel = get_state_model(
                don.get(index), stt.get(index), descr='signal')
        donBlcModel = get_state_model(
                don.get(index), stt.get(index), descr='acceptor_bleached')
        accModel = get_state_model(
                acc.get(index), stt.get(index), descr='signal')
        fretModel = get_state_model(
                fret.get(index), stt.get(index), descr='signal')

        fig.donModelOverlay.draw(donSigModel, color='black')
        fig.donModelOverlay.draw(donBlcModel, color='black')
        fig.accModelOverlay.draw(accModel, color='black')
        fig.fretModelOverlay.draw(fretModel, color='red')

        if 'clusters' in descr:
            clusters = traces[descr.index('clusters')]
            clusModel = clusters.get(index)
            fig.clusModelOverlay.draw(clusModel, color='blue')
        
        if self.get_item('model').get()==False:
            fig.donModelOverlay.toggle(False)
            fig.accModelOverlay.toggle(False)
            fig.fretModelOverlay.toggle(False)
            fig.clusModelOverlay.toggle(False)


        # Bleach event overlay
        accBleachCrit = stt.data[index,:] == stt.codes['acceptor_bleach_event']
        accBleachEvent = acc.time[accBleachCrit]
        donBleachCrit = stt.data[index,:] == stt.codes['donor_bleach_event']
        donBleachEvent = don.time[donBleachCrit]

        dt = np.mean(np.diff(don.time))
        h = []
        if donSigModel.data.size > 0:
            h.append(np.max(donSigModel.data))
        if donBlcModel.data.size > 0:
            h.append(np.max(donBlcModel.data))
        if accModel.data.size > 0:
            h.append(np.max(accModel.data))
        
        if len(h) > 0:
            h = 1.2*np.max(h)

        if accBleachEvent.size > 0:
            accBleachTime = np.mean(accBleachEvent)
            fig.accBleachOverlay.draw(accBleachTime, h, width=3*dt, facecolor='red')

        if donBleachEvent.size > 0:
            donBleachTime = np.mean(donBleachEvent)
            fig.donBleachOverlay.draw(donBleachTime, h, width=3*dt, facecolor='green')

        fig.show()

    def clear(self, *args, **kwargs):
        fig = self.tracePanel.figureHandle
        if self.tracePanel.has_window():
            fig.donModelOverlay.clear()
            fig.accModelOverlay.clear()
            fig.fretModelOverlay.clear()
            fig.clusModelOverlay.clear()

            fig.accBleachOverlay.clear()
            fig.donBleachOverlay.clear()


class ModelOverlay():
    """ Overlay HMM model on trace. """

    def __init__(self, ax=None):
        self.ax = ax
        self.overlay = []

    def draw(self, data, *args, **kwargs):
        if data.time.shape[0] == data.data.shape[1]:
            dData = np.transpose(data.data)
        else:
            dData = data.data
        p = self.ax.plot(data.time, dData, *args, **kwargs)
        self.overlay.append(p[0])

    def clear(self):
        try:
            [p.remove() for p in self.overlay]
        except:
            pass
        self.overlay = []
    
    """
    def update(self, index, *args, **kwargs):
        if not index==None:
            if index==-1:
                [p.set(*args, **kwargs) for p in self.overlay]
            else:
                self.overlay[index].set(*args, **kwargs)
    """

    def toggle(self, visible):
        if type(visible) == bool:
            visible = visible*np.ones(len(self.overlay))
        [self.overlay[m].set_visible(visible[m]) for m in range(len(self.overlay))]


class BleachZoneOverlay():
    """ Overlay indicator for bleaching events. """
    
    def __init__(self, ax=None):
        self.ax = ax
        self.overlay = []

    def draw(self, *args, **kwargs):
        #p = self.ax.plot(data.time, data.data, *args, **kwargs)
        p = self.ax.bar(*args, edgecolor=None, alpha=0.5, align='center', **kwargs)
        d = drag(p[0], 
                allow_resize=False,
                fixed_aspect_ratio=False,
                lock_y=True)
        d.connect()

        self.overlay.append(d)

    def clear(self):
        try:
            [p.rect.remove() for p in self.overlay]
        except:
            pass
        self.overlay = []

    def toggle(self, visible):
        if type(visible) == bool:
            visible = visible*np.ones(len(self.overlay))
        [self.overlay[m].rect.set_visible(visible[m]) 
                for m in range(len(self.overlay))]

