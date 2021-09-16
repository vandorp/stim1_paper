""" Timeseries filters. """

import numpy as np
from math import factorial


# Filter definitions ----------------------------------------------------------
def moving_average(data, winSize):
    if np.prod(data.shape) == data.size:
        data = np.reshape(data, data.size)

    weights = np.repeat(1.0, winSize)/winSize
    return np.convolve(data, weights, 'same')


def moving_rms(data, winSize):
    data = data - moving_average(data, winSize)
    return np.sqrt(moving_average(data*data, winSize))

#*#
def savitzky_golay(data, winSize, order, deriv=0, rate=1):
    order_range = range(order+1)
    half_window = (winSize -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, 
        half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = data[0] - np.abs(data[1:half_window+1][::-1] - data[0])
    lastvals = data[-1] + np.abs(data[-half_window-1:-1][::-1] - data[-1])
    data = np.concatenate((firstvals, data, lastvals))
    try:
        return np.convolve( m[::-1], data, mode='valid')
    except:
        return data


# Classes ---------------------------------------------------------------------
#*#
class FIRFilter(object):
    """ Base class for FIR filters. """

    def __init__(self, func, args):
        self.func = func
        self.args = args
        
    def run(self, ts):
        data = ts.data
        nrTraces = ts.get_nr_traces()
        if nrTraces > 1:
            for n in range(nrTraces):
                data[n,:] = self.func(data[n,:], *self.args)
        else:
            data = self.func(data, *self.args)
        ts.data = data

#*#
class MovingWindow(FIRFilter):
    """ General moving window. """

    def __init__(self, func, win, *args):
        
        try:
            win = np.abs(np.int(win))
        except ValueError:
            raise ValueError("Window size has to be of type int")
        if win % 2 != 1 or win < 1:
            raise TypeError("Window size size must be a positive odd number")

        FIRFilter.__init__(self, func, [win]+list(args))

    def get_window(self):
        return self.args[0]


class MovingAverage(MovingWindow):
    """ Average value in moving window. """
    
    def __init__(self, win=1):
        MovingWindow.__init__(self, moving_average, win)


class MovingRMS(MovingWindow):
    """ Average RMS in moving window. """

    def __init__(self, win=1):
        MovingWindow.__init__(self, moving_rms, win)

#*#
class SavGol(MovingWindow):
    """ Polynomial fit in moving window. """
    
    def __init__(self, windowSize, order):
        win = windowSize
        try:
            order = np.abs(np.int(order))
        except ValueError:
            raise ValueError("Order has to be of type int")

        if win < order + 2:
            raise TypeError("Window size is too small for the polynomials \
                    order")

        MovingWindow.__init__(self, savitzky_golay, win, order)

