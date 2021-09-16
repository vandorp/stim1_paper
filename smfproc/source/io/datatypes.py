""" Datatype definitions. """

import numpy as np


""" ===========================================================================
= MOVIES
===============================================================================
"""

class ImageCoordinates():
    """ Storing 2D coordinates of recorded molecules in a movie.

    Each molecule receives (row, column) coordinates according to its location
    in a movie frame.
    """

    def __init__(self, crds=[], imageShape=(0,0)):
        self.imageShape = imageShape
        
        if isinstance(crds, list):
            self.crds = crds
        
        elif isinstance(crds, ImageCoordinates):
            self.crds = list(crds.crds)
            self.imageShape = tuple(crds.imageShape)
            #self = ImageCoordinates(crds.crds, crds.imageShape)

        else:
            print('Invalid data type in Coordinates')
            self.crds = []
            self.imageShape = (0,0)

    def select(self, crit):
        crds = []
        for i in range(len(self.crds)):
            if crit[i]:
                crds.append(self.crds[i])
        return ImageCoordinates(crds, self.imageShape)

    def get_crd(self, index):
        return ImageCoordinates([self.crds[index]], self.imageShape)
        #return self.select(index)

    def get(self):
        return self.crds

    def set(self, crds=[]):
        self.crds = crds

    def get_row_col(self):
        m = lambda n: np.array([c[n] for c in self.crds], dtype=int)
        return [m(0), m(1)]

    def set_row_col(self, row, col):
        self.crds = [(row[n], col[n]) for n in range(len(row))]

    def get_linear(self):
        pass

    def set_linear(self):
        pass

    def get_distances(self):
        nrCoords = len(self.crds)
        M = np.zeros((nrCoords, nrCoords))
        for n in range(nrCoords):
            for m in range(nrCoords):
                M[n,m] = np.sqrt(np.sum((np.array(self.crds[n]) - np.array(self.crds[m]))**2))
        return M


class Movie():
    """ Movie data.

    Movie frames are stored as 3D numpy arrays.
    """
    
    def __init__(self, data, header):

        #self.time = time
        self.header = header
        
        if isinstance(data, np.ndarray):
            # Data must be organized as (rows, cols, frames)
            self.frames = data
                   
        elif isinstance(data, SmfMovie):
            self = Movie(data.frames, data.dt, data.header)

        else:
            print('Invalid data type in SmfMovie')
            self.frames = np.array([])

    def width(self):
        return self.get_frame(0).width()

    def height(self):
        return self.get_frame(0).height()

    def shape(self):
        return self.get_frame(0).shape()

    def nr_frames(self):
        return self.frames.shape[2]

    def split(self):
        return [Movie(s, self.header) for s in np.hsplit(self.frames, 2)]

    def get_average_frame(self, indices=-1):
        if indices == -1:
            frames = self.frames
        else:
            frames = self.frames[:,:,indices]
        return Frame(np.mean(frames, axis=2))

    def apply(self, func, *args, **kwargs):
        return func(self.frames, *args, **kwargs)

    def get_frame(self, index):
        return Frame(self.frames[:,:,index])

    def get(self, index=-1):
        if index==-1:
            return self.frames
        else:
            return self.frames[:,:,index]
    
    def get_header_str(self):
        hstr = ['%s: %s'%(key, self.header[key]) for key in self.header]
        return '\n'.join(hstr)


class Frame():
    """ Movie frame.
    
    Frames are stored as 2D numpy arrays.
    """

    def __init__(self, data):

        if isinstance(data, np.ndarray):
            self.pixels = np.array(data, dtype=float)

        elif isinstance(data, SmfImage):
            self = Frame(data.pixels) # This may not work as expected!

        else:
            print('Invalid data type in SmfImage')
            self.pixels = np.array([])

    def width(self):
        return self.pixels.shape[1]

    def height(self):
        return self.pixels.shape[0]

    def shape(self):
        return self.pixels.shape

    #def split(self):
    #    return [SmfImage(s) for s in np.hsplit(self.pixels, 2)]

    def apply(self, func, *args, **kwargs):
        return func(self.pixels, *args, **kwargs)
    
    def get(self):
        return self.pixels

    def set(self, pixels):
        self.pixels = pixels



""" ===========================================================================
= TIMESERIES
===============================================================================
"""

class TimeSeries(object):
    """ Storing and manipulating timeseries data.
    """
    
    def __init__(self, data, time=np.array([])):
        
        #if isinstance(data, TimeSeries):
        try:
            self.data = np.array(data.data)
            self.time = np.array(data.time)
        except:
            self.data = np.array(data)

            #if time==None: # NB: 'None' is not a good way to do this 
            if len(time)==0:
                if len(self.data.shape) > 0:
                    self.time = np.arange(self.data.shape[-1])
                else:
                    self.time = np.array([])
            else:
                self.time = np.array(time)
        
        if len(self.data.shape) == 1:
            self.data = np.reshape(self.data, [1, self.data.size])
        elif len(self.data.shape) > 1:
            if self.data.shape[0] == self.time.size:
                self.data = np.transpose(self.data)
        
    def get_nr_traces(self):
        nrRows = self.data.shape[0]
        if len(self.data.shape)==1:
            if self.data.shape[0]==0:
                nrRows = 0
            else:
                nrRows = 1
        
        if self.data.size == np.max(self.data.shape):
            nrRows = 1
        
        if self.data.size == 0:
            nrRows = 0

        return nrRows

    def get(self, index):
        if self.get_nr_traces() == 1:
            return TimeSeries(self.data, self.time)
        elif self.get_nr_traces() == 0:
            return TimeSeries([])
        else:
            data = self.data[index,:]
            #data = np.reshape(data, [1,data.size])
            return TimeSeries(data, self.time)

    def select(self, crit):
        if len(crit)==0:
            return TimeSeries(np.array([]))
        else:
            data = self.data[crit==1,:]
            #data = np.reshape(data, [1,data.size])
            return TimeSeries(data, self.time)

    def slice(self, crit):
        crit = np.reshape(crit, crit.size)
        time = np.reshape(self.time, self.time.size)
        if self.get_nr_traces() < 2:
            data = np.reshape(self.data, self.data.size)
            return TimeSeries(data[crit], time[crit])
        else:
            return TimeSeries(self.data[:, crit], time[crit])

    def stack(self, ts):
        data = np.vstack([self.data] + [d.data for d in ts])
        time = self.time
        return TimeSeries(data, time)

    def histogram(self, *args, **kwargs):
        cnts = []
        bins = []
        for n in range(self.get_nr_traces()):
            cnt, binEdges = np.histogram(self.data[n,:], *args, **kwargs)
            cnts.append(cnt)
            bins = (binEdges[:-1] + binEdges[1:])/2
        return Histogram(cnts, bins)

    def copy(self):
        return TimeSeries(self)

    def plot(self, ax, *args, **kwargs):
        if self.get_nr_traces() == 1:
            self.time = np.reshape(self.time, self.time.size)
            self.data = np.reshape(self.data, self.data.size)

        #if self.data.shape.index(self.time.size) == 1:
        #    self.data = np.transpose(self.data)
        data = self.data
        #data = data[np.isnan(data)==False]
        #time = self.time[np.isnan(data)==False]
        data[np.isnan(data)] = 0

        return ax.plot(self.time, data, *args, **kwargs)


class State(TimeSeries):
    """ Derived timeseries class for representing idealized processes, 
    generated e.g. by an HMM filter.
    """

    codes = {
            'misc': 0,
            'acceptor_bleach_event': 2,
            'donor_bleach_event': 3,
            'signal': 10,
            'acceptor_bleached': 20,
            'donor_blinked': 29,
            'donor_bleached': 30,
            'acceptor_present': 40}

    def __init__(self, data, time=np.array([])):
        TimeSeries.__init__(self, data, time)
        self.data = np.array(self.data, dtype=int)

    def select(self, crit):
        return State(TimeSeries.select(self, crit))

    def has(self, descr):
        code = self.codes[descr]
        if self.get_nr_traces() > 1:
            return np.sum(self.data == code, axis=1) > 0
        else:
            return np.sum(self.data == code) > 0

    def set_transition(self, index, n, descr, zone=3):
        code = self.codes[descr]
        eventCode = int(code/10)

        if zone > 0:
            w = int((zone-1)/2)+1
            self.data[index, n-w+1:n+w] = eventCode
        else:
            w = 0
        self.data[index, n+w:] = code
    
    def get_codes(self, descr=None):
        if descr==None:
            return sorted(list(set(self.data)))
        else:
            crit = self.get_all_crit(descr)
            data = self.slice(crit).data
            data = np.reshape(data, data.size)
            return sorted(list(set(data)))

    def get(self, index):
        return State(TimeSeries.get(self, index))

    def get_crit(self, index=0, descr=None):

        data = TimeSeries.get(self, index)
        
        if descr==None:
            codes = sorted(list(set(self.data)))
        else:
            code = self.codes[descr]

            crit = (data.data >= code) & (data.data < code+10)
            stateData = data.slice(crit)
            
            stData = np.reshape(stateData.data, stateData.data.size)
            codes = sorted(list(set(stData)))
        allData = []
        for c in codes:
            allData.append(data.data==c)

        return allData

    def get_all_crit(self, descr):
        code = self.codes[descr]
        return (self.data >= code) & (self.data < code+10)

    def set_crit(self, index, crit, descr):
        code = self.codes[descr]
        
        """
        if crit.size in crit.shape:
            crit = np.reshape(crit, (1, crit.size))
            nrRows = 1
        else:
            nrRows = crit.shape[0]

        for n in range(nrRows):
            self.data[index, crit[n,:]==1] = code + n
        """

        if not (type(crit)==list):
            crit = [crit]
        
        for n in range(len(crit)):
            self.data[index, crit[n]==1] = code + n


    def set_all_crit(self, crit, descr):
        self.data[crit==1] = self.codes[descr]



class SmfSpecs(object):
    """ Storing TimeSeries metadata.
    """

    def __init__(self, data=np.array([]), descr=list()):
        #if isinstance(data, SmfSpecs):
        try:
            self.descr = data.descr
            self.data = data.data
        except:
            #else:
            # data must be column-wise
            self.descr = list(descr)
            data = np.array(data)
            if len(data.shape) == 1:
                data = np.reshape(data, [1, data.size])
            self.data = data

    def get(self, descr=None):
        if descr==None:
            data = self.data
        else:
            if descr in self.descr:
                index = self.descr.index(descr)
                data = self.data[:,index]
            else:
                data = np.array([])
        
        return data

    def get_nr_rows(self):
        return self.data.shape[0]

    def get_nr_columns(self):
        return self.data.shape[1]

    def select(self, crit):
        data = self.data[crit,:]
        return SmfSpecs(data, self.descr)

    def set(self, data, descr):
        data = np.array(data)
        if descr in self.descr:
            index = self.descr.index(descr)
            self.data[:,index] = data
        else:
            self.descr.append(descr)
            data = np.reshape(data,(data.size,1))
            if (self.get_nr_columns() > 0) & (self.get_nr_rows() > 0):
                self.data = np.hstack((self.data, data))
            else:
                self.data = data

    def merge(self, specs):
        nrColumns = specs.get_nr_columns()
        for n in range(nrColumns):
            self.set(specs.data[:,n], specs.descr[n])

    def append_old(self, specs):
        if self.data.size > 0:
            if cmp(self.descr, specs.descr)==0:
                self.data = np.vstack((self.data, specs.data))
            else:
                print('Error. Column descriptions must match.')
        else:
            self.__init__(specs)

    def append(self, specs):
        if self.data.size > 0:
            descrInOldButNotInNew = list(set(self.descr) - set(specs.descr))
            descrInNewButNotInOld = list(set(specs.descr) - set(self.descr))
            descrOverlap = list(set(self.descr) & set(specs.descr))

            #print descrInOldButNotInNew
            #print descrInNewButNotInOld
            #print descrOverlap

            nrOldData = self.data.shape[0]
            nrNewData = specs.data.shape[0]

            #cols = np.array((0,nrOldData+nrNewData))
            cols = np.hstack((
                self.data[:,self.descr.index(descrOverlap[0])],
                specs.data[:,specs.descr.index(descrOverlap[0])]))

            for d in descrOverlap[1:]:
                col = np.hstack((
                        self.data[:,self.descr.index(d)],
                        specs.data[:,specs.descr.index(d)]))
                cols = np.vstack((cols, col))

            for d in descrInOldButNotInNew:
                col = np.hstack((
                        self.data[:,self.descr.index(d)],
                        np.zeros(specs.data.shape[0])))
                cols = np.vstack((cols, col))

            for d in descrInNewButNotInOld:
                col = np.hstack((
                        np.zeros(self.data.shape[0]),
                        specs.data[:,specs.descr.index(d)])),
                cols = np.vstack((cols, col))

            #self.data = np.transpose(np.array(cols))
            self.data = np.transpose(cols)
            self.descr = descrOverlap + descrInOldButNotInNew + descrInNewButNotInOld

        else:
            self.__init__(specs)


class SmfSelect(SmfSpecs):
    """ Specialized metadata class for specifying selected traces.
    """

    def __init__(self, data=np.array([]), descr=list()):
        if isinstance(data, SmfSpecs):
            SmfSpecs.__init__(self, data.data, data.descr)
        else:
            SmfSpecs.__init__(self, data, descr)
        self.data = np.array(self.data) > 0

    def set(self, data, descr):
        data = np.array(data) > 0
        SmfSpecs.set(self, data, descr)

    def get(self, descr=None, func=None):
        if func==None:
            return SmfSpecs.get(self, descr)
        else:
            if len(descr) < 2:
                return SmfSpecs.get(self, descr[0])
            else:
                A = SmfSpecs.get(self, descr[0])
                for n in range(1,len(descr)):
                    A = func(A, SmfSpecs.get(self, descr[n]))
                return A


