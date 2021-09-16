import matplotlib.pyplot as plt
from pylab import get_current_fig_manager


class FigureBase(object):
    """ Base class for output figures. """

    def __init__(self, *args, **kwargs):
        self.figure = plt.figure(*args, **kwargs)

    def set_position(self, crds):
        plt.figure(self.figure.number)
        man = get_current_fig_manager()
        man.window.wm_geometry("+%d+%d"%(crds[0],crds[1]))

    def set_size(self, figSize):
        self.figure.set_size_inches(figSize)

    def set_title(self, title):
        if not (self.figure==None):
            self.figure.canvas.set_window_title(title)

    def show(self):
        self.figure.show()

    def close(self):
        plt.close(self.figure)

