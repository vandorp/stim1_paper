""" Base classes for GUI components. """

import os
import numpy as np
import matplotlib.pyplot as plt

from . import Tk


class VarEntry(object):
    """ For entry of variables. """

    def __init__(self, master, label):
        self.var = Tk.StringVar()
        self.frame = Tk.Frame(master)
        self.label = Tk.Label(master=self.frame, text=label)
        self.entry = None
    
    def pack(self, *args, **kwargs):
        self.label.pack(side=Tk.LEFT, fill=Tk.BOTH)
        self.entry.pack(side=Tk.RIGHT, fill=Tk.BOTH)
        self.frame.pack(*args, **kwargs)

    def config(self, *args, **kwargs):
        self.entry.config(*args, **kwargs)

    def get(self):
        return self.var.get()

    def set(self, value):
        self.var.set(value)


class BoolEntry(VarEntry):
    """ For entry of boolean variables. """

    def __init__(self, master, label='', *args, **kwargs):
        VarEntry.__init__(self, master, label)
        self.entry = Tk.Checkbutton(master=self.frame, 
                variable=self.var, *args, **kwargs)

    def select(self):
        return self.entry.select()

    def deselect(self):
        return self.entry.deselect()

    def get(self):
        return VarEntry.get(self)=='1'


class NumEntry(VarEntry):
    """ For entry of numeric values. """

    def __init__(self, master, label=''):
        VarEntry.__init__(self, master, label)
        self.entry = Tk.Entry(master=self.frame, textvariable=self.var, width=4)
        self.var.set(0)

    def get(self):
        strValue = VarEntry.get(self)
        if strValue.find('.') == -1:
            return int(strValue)
        else:
            return float(strValue)

    def set(self, value):
        self.var.set(value)



class UpdateCallback(object):
    """ Base class for handling GUI callbacks. """

    def __init__(self):
        self.callbacks = []

    def add_update_callback(self, func_def, loc=0):
        #self.callbacks.append(func_def)
        index = len(self.callbacks) - loc
        self.callbacks.insert(index, func_def)

    def callback(self, val=None, *args):
        [func(val) for func in self.callbacks]


class TkItems(object):
    """ Base class for storing references to GUI panels. """

    def __init__(self):
        self.items = []
        self.descr = []
        self.vars = []

    def add_item(self, item, descr, var=None):
        #self.items[descr] = [item, var]
        self.items.append(item)
        self.descr.append(descr)
        self.vars.append(var)

    def get_item(self, descr):
        return self.items[self.descr.index(descr)]

    def get_var(self, descr):
        return self.vars[self.descr.index(descr)]
    

class WindowPanel(Tk.Toplevel, TkItems):
    """ Base class for top level window. """

    def __init__(self, master, descr=''):
        Tk.Toplevel.__init__(self, master)
        TkItems.__init__(self)
        
        master.add_item(self, descr)
        self.master = master

        self.title(descr)
        self.hidden = False
    
    def hide(self):
        self.withdraw()
        self.hidden = True

    def show(self):
        self.deiconify()
        self.hidden = False

    def pack(self):
        #Panel.pack(self, side=Tk.LEFT, fill=Tk.BOTH)
        self.config(relief='groove', borderwidth=1)
        [item.pack() for item in self.items]

    def set_position(self, crds):
        self.geometry('+%d+%d'%tuple(crds))

    def _quit(self):
        pass


class Panel(Tk.Frame, TkItems, UpdateCallback):
    """ Base class for GUI panels. """
    
    def __init__(self, master, descr):
        Tk.Frame.__init__(self, master=master)
        TkItems.__init__(self)
        UpdateCallback.__init__(self)

        master.add_item(self, descr) # This seems superfluous...
        self.master = master

    def _quit(self):
        for item in self.items:
            try:
                item._quit()
            except Exception as e:
                #print e.message
                pass


class SubPanel(Panel):
    """ GUI subpanel. """

    def pack(self):
        Panel.pack(self, side=Tk.TOP, fill=Tk.BOTH)
        self.config(relief='groove', borderwidth=1)
        [item.pack(side=Tk.TOP, fill=Tk.BOTH) for item in self.items]

class EntryPanel(SubPanel):
    """ GUI subpanel for value entries. """

    def set(self, descr, value):
        self.get_item(descr).set(value)

    def set_group(self, values=dict()):
        for key in values:
            self.set(key, values[key])

    def get(self, descr):
        return self.get_item(descr).get()



class TopPanel(Panel):
    """ Main panel. """

    def __init__(self, master, descr):
        Panel.__init__(self, master, descr)
        
        label = Tk.Label(master=self, text=descr, font=("Helvetica", 12))
        self.add_item(label, 'label')

    def pack(self):
        Panel.pack(self, side=Tk.LEFT, fill=Tk.BOTH)
        self.config(relief='groove', borderwidth=1)
        [item.pack() for item in self.items]
    

class FigurePanel(SubPanel):
    """ Panel for controlling output figures. """

    def __init__(self, master, descr):
        SubPanel.__init__(self, master, descr)

        self.data = []
        self.figureHandle = None
        self.figurePosition = [0,0]
        self.figureSize = [5,5]

    def has_data(self):
        return len(self.data) > 0

    def has_window(self):
        return not (self.figureHandle == None)
    
    def open_window(self, fclass):
        self.figureHandle = fclass(figsize=self.figureSize)
        #self.figureHandle.set_size(self.figureSize)
        self.figureHandle.set_position(self.figurePosition)

    def close_window(self):
        if self.has_window():
            self.figureHandle.close()
            self.figureHandle = None


def get_monitor_size():
    """ Returns pixel dimensions of the monitor """

    width = Tk.Tk().winfo_screenwidth()
    height = Tk.Tk().winfo_screenheight()
    return width, height


class ScrollBarGUI(Tk.Frame, UpdateCallback):
    """ For making a GUI scrollbar object. """

    def __init__(self, master):
        Tk.Frame.__init__(self, master=master)
        UpdateCallback.__init__(self)

        self.scrollbar = Tk.Scale(
                master=self, 
                command=self.update,
                orient=Tk.HORIZONTAL, 
                showvalue=False)

        self.max_nr = 0
            
    def pack(self, *args, **kwargs):
        self.scrollbar.pack(fill=Tk.X)
        Tk.Frame.pack(self, *args, **kwargs)

    def config(self, nr):
        self.max_nr = nr
        self.scrollbar.config(
                from_=0,
                to=nr-1,
                resolution=1)
        self.reset()

    def update(self, *args):
        val = int(args[0])
        self.callback(val)

    def reset(self):
        self.scrollbar.set(0)

    def get(self):
        return self.scrollbar.get()

    def set(self, *args):
        self.scrollbar.set(*args)

    def set_last(self):
        self.set(self.max_nr-1)

    def set_first(self):
        self.set(0)

    def key_prev(self, event):
        cval = self.scrollbar.get()
        nval = cval - 1
        if cval > 0:
            self.scrollbar.set(nval)

    def key_next(self, event):
        nr = self.max_nr
        cval = self.scrollbar.get()
        nval = cval + 1
        if cval < (nr-1):
            self.scrollbar.set(nval)


class ListBox(UpdateCallback):
    """ For making a GUI listbox. """

    def __init__(self, master, **kwargs):
        UpdateCallback.__init__(self)

        self.listBox = Tk.Listbox(
                master=master,
                exportselection=0,
                background='white',
                selectbackground='gray75', 
                **kwargs)
        self.listBox.bind("<<ListboxSelect>>", self.update)

        self.execute_at_end = []
        self.execute_at_start = []

        self.subset = np.array([])
        #self.currentIndex = -1

        #def get_subset(self):
        #return self.subset

        self.content = []
        self.display = []
        
    def next(self, *args, **kwargs):
        current = self.listBox.curselection()[0]
        self.listBox.selection_clear(current)
        self.listBox.selection_set(current+1)
        self.update()
    
    def prev(self, *args, **kwargs):
        current = self.listBox.curselection()[0]
        self.listBox.selection_clear(current)
        self.listBox.selection_set(current-1)
        self.update()

    def set(self, content=[], display=[]):
        self.clear()

        self.content = content

        if len(display)==0:
            self.display = content
        else:
            self.display = display

        self.display = ['%s'%(self.display[n]) 
                for n in range(len(self.display))]
        
        #self.show()
    
    def get(self, index=None):
        #return self.currentIndex
        if self.size() > 0:
            if index==None:
                currentIndex = self.listBox.curselection()[0]
            else:
                currentIndex=index
            currentContent = self.listBox.get(currentIndex)
            #index = strcmp(self.display, currentContent)
            index = np.array([display==currentContent for display in self.display])==1
            return np.array(self.content)[index][0]
        else:
            return None

    def get_true_index_from_subset_index(self, subsetIndex):
        displayContent = self.listBox.get(subsetIndex)
        indexCrit = np.array([display==displayContent 
            for display in self.display])
        return np.where(indexCrit==1)
    
    def pack(self, *args, **kwargs):
        self.listBox.pack(*args, **kwargs)

    def clear(self):
        display = self.listBox
        display.delete(0, Tk.END)
        self.currentIndex = -1
 
        self.listBox.config(
                background='white',
                selectbackground='gray75')

    def show(self, subset=np.array([])):
        subset = np.array(subset)

        if subset.size==0:
            subset = np.ones(len(self.content))==1

        if len(subset.shape)==0:
            subset = np.reshape(subset, subset.size)
        
        self.subset = subset==1
        #self.displayed = subset
        
        #print self.display, self.subset
        values = np.array(self.display)[self.subset]
        self.clear()
        lb = self.listBox
        #display.delete(0, Tk.END)
        [lb.insert(Tk.END, value) for value in values]
        
        if len(values) > 0:
            self.select(0)
            #self.currentIndex = 0

    def size(self):
        return self.listBox.size()

    def update(self, *args, **kwargs):
        #if not (self.currentIndex < 0):
        #    oldValue = self.get(self.currentIndex)
        #else:
        #    oldValue = None
        value = self.get()

        if self.size() > 0:
            #index = self.listBox.curselection()[0]

            #if not (value == oldValue):
            #    self.currentIndex = index
            self.callback(value)

    def select(self, index):
        self.listBox.selection_set(index)


class SuperListBox(ListBox):
    """ For making a listbox with additional selection options. """

    def __init__(self, master):
        ListBox.__init__(self, master)
        
        self.mode = np.array([]) # To vary e.g. selection color
        self.superSelected = np.array([])
        self.superSelectedCallback = lambda:None
        
        self.listBox.bind('<Double-Button-1>', self.super_select_toggle)
        self.listBox.bind('<space>', self.super_select_toggle)

    def take_focus(self):
        self.listBox.focus_force()

    def set(self, content=[], display=[]):
        ListBox.set(self, content, display)
        self.superSelected = np.zeros(len(self.display))

    def show(self, subset=np.array([])):
        ListBox.show(self, subset)
        
        """
        sup = self.superSelected

        if sup.size == 0:
            superDisplay = np.array([])
        elif sup.size == 1:
            if self.subset == 0:
                superDisplay = np.array([])
            elif self.subset == 1:
                superDisplay = np.array(sup)
        elif sup.size > 1:
        """
        
        A = np.array(self.superSelected)
        B = np.reshape(A, A.size)

        superDisplay = B[self.subset]
        
        for d in range(len(superDisplay)):
            s = superDisplay[d]
            #if ((s==1) | (s==True)):
            index = self.get_true_index_from_subset_index(d)
            if s:
                mode = self.mode[index][0]
                self.super_select_index(d, mode)

    def clear(self):
        ListBox.clear(self)
        #self.superSelected = np.array([])
    
    def super_select_toggle(self, *args, **kwargs):
        displayIndex = self.listBox.curselection()[0]
        #displayContent = self.listBox.get(displayIndex)
        #indexCrit = np.array([display==displayContent 
        #    for display in self.display])
        #index = np.where(indexCrit==1)
        index = self.get_true_index_from_subset_index(displayIndex)
        
        #np.array(self.content)[index][0]

        if self.item_is_super_selected(index):
            self.superSelected[index] = False
            self.super_deselect_index(displayIndex)
        else:
            self.superSelected[index] = True
            mode = self.mode[index][0]
            self.super_select_index(displayIndex, mode)

        #self.set_super_selected(indices)
        self.superSelectedCallback()

    def super_select_index(self, index, mode=False):
        if mode:
            self.listBox.itemconfig(index, 
                    background='red', 
                    selectbackground='orange')
        else:
            self.listBox.itemconfig(index, 
                    background='turquoise', 
                    selectbackground='steelblue3')

    def super_deselect_index(self, index):
        self.listBox.itemconfig(index, 
                background='white', 
                selectbackground='gray75')

    def item_is_super_selected(self, index):
        #bg = self.listBox.itemcget(index, 'background')
        #noBg = (len(bg) == 0)
        #whiteBg = (bg == 'white')
        #return not (noBg | whiteBg)
        return self.superSelected[index]

    def get_super_selected(self):
        #return [self.is_super_selected(index) 
        #        for index in range(self.listBox.size())]
        return self.superSelected

    def set_super_selected(self, crit):
        #[self.super_select(index) for index in indices]
        self.superSelected = np.array(crit)==1


class OptionMenu(Tk.OptionMenu):
    """ Makes a dropdown menu. """
    
    def __init__(self, *args, **kw):
        self._command = kw.get("command")
        Tk.OptionMenu.__init__(self, *args, **kw)
        self.var = args[1]
    
    def add_option(self, label, variable):
        self["menu"].add_command(label=label,
                command=Tk._setit(variable, label, self._command))

    def clear(self):
        self["menu"].delete(0,Tk.END)

    def size(self):
        return len(self.var.get())

"""
import re
class HoverInfo(Tk.Menu):
    
    def __init__(self, parent, text, command=None):
        self._com = command
        Tk.Menu.__init__(self,parent, tearoff=0)
        if not isinstance(text, str):
            raise TypeError(
                    'Trying to initialise a Hover Menu with a non string type: ' 
                    + text.__class__.__name__)
        toktext=re.split('\n', text)
        for t in toktext:
            self.add_command(label = t)
        self._displayed=False
        self.master.bind("<Enter>",self.Display )
        self.master.bind("<Leave>",self.Remove )

    def __del__(self):
        self.master.unbind("<Enter>")
        self.master.unbind("<Leave>")

    def Display(self,event):
        if not self._displayed:
            self._displayed=True
            self.post(event.x_root, event.y_root)

        if self._com != None:
            self.master.unbind_all("<Return>")
            self.master.bind_all("<Return>", self.Click)

    def Remove(self, event):
        if self._displayed:
            self._displayed=False
            self.unpost()

        if self._com != None:
            self.unbind_all("<Return>")

    def Click(self, event):
        self._com()
"""
