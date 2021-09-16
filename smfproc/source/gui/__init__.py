from matplotlib import use as useBackend
useBackend('TkAgg')

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    import tkFileDialog
else:
    import tkinter as Tk
    from tkinter import filedialog as tkFileDialog
    #tkFileDialog = Tk.filedialog #import pysmf_v3.core as core
