""" Data file definitions. """



# System packages
import os, stat, sys
import numpy as np
import warnings
import datetime
import struct

# Program modules
from smfproc.source.io.datatypes import Movie, TimeSeries, SmfSpecs, SmfSelect


""" ===========================================================================
= MOVIES
===============================================================================

Supported movie files are stacked tiff (.tiff) as exported by uManager, and an 
old custom format from TJ Ha lab (.pma). Movie files are read-only.

"""

class MovieFile():
    """ Generic interface for movie data.

    Determines the movie type based on the file extension, and then casts
    the corresponding class onto this generic class, so that all movie data
    can be accessed in the same way by the program.
    """

    def __init__(self, filePath):
        
        movieFile = os.path.sep.join(filePath.split(os.path.sep)[-2:])
        print('Initializing movie %s'%movieFile)
        
        self.filePath = filePath
        #self.meta = {
        #            'summary':{},
        #            'frames':{},
        #            'timeseries':lambda : None}

        self.header = {}

        movieExt = filePath.split('.')[-1]
        if movieExt == 'tif':
            movie = MMstackedTIFF(filePath)

            meta = movie.get_metadata()
            summary, frames = meta.read()
            #self.meta['summary'] = summary
            #self.meta['frames'] = frames
            #self.meta['timeseries'] = meta.get_timeseries
            
            self.header['dt'] = \
                    float(movie.log['Andor-ActualInterval-ms'])/1000.0
            self.header['imageWidth'] = movie.header['ImageWidth']
            self.header['imageHeight'] = movie.header['ImageHeight']
            self.header['nrFrames'] = movie.NrFrames

            self.read = movie.read
            self.get_timeseries = meta.get_timeseries

            
        elif movieExt == 'pma':
            movie = PmaFile(filePath)
            self.header['imageWidth'] = movie.header["WIDTH"]
            self.header['imageHeight'] = movie.header["HEIGHT"]
            self.header['nrFrames'] = movie.header["NR_FRAMES"]
            self.header['dt'] = 0.1
            self.read = movie.read
            self.get_timeseries = lambda *args: np.array([])

     


class MMstackedTIFF():
    """
    TIFF file starts with an 8-byte file header containing:
        - Byte order ('II' or 'MM')
        - A signature (The number 42)
        - Location of the first IFD

    Every image in the file is preceded by an Image File Directory (IFD). The
    IFD contains info about image properties, and is directly followed by the
    image data. Each IFD also points to the location of the next IFD.
    """

    encoding = 'utf-8'

    IFDtags = {
            256: 'ImageWidth', 
            257: 'ImageHeight',
            258: 'BitsPerSample',
            51123: 'LogLoc'}
    
    metaTag = '_metadata.txt'

    def __init__(self, FilePath):
        self.FilePath = FilePath

        self.NrFrames = None
        self.DataLocs = []
        self.header = self._parse_headers()
        self.log = self._parse_log()
        #print self.log
        #self._parse_metadata()
        #self.read_data()
        
        ijmeta = self.get_metadata()
        
        #summary,_ = self.get_metadata()
        #self.ImageHeight = int(summary['Height'])
        #self.ImageWidth = int(summary['Width'])
        self.NrFrames = int(ijmeta.summary['Frames'])

    def get_metadata(self):
        metaFile = self.FilePath.split('.ome')[0] + self.metaTag
        #summary = {}
        #frames = {}
        if os.path.isfile(metaFile):
            return IJMetaFile(metaFile)
            #summary, frames = M.read()
        else:
            return None
        #return summary, frames
        
    def _parse_headers(self):
        fid = open(self.FilePath, 'rb')

        """ File header """
        # 'II' = little-endian, 'MM' = big-endian
        self.ByteOrder = struct.unpack('2c', fid.read(2))
        # TIFF sanity check. This should of course yield the number 42
        TheNumber42 = struct.unpack('H', fid.read(2))[0]
        if not TheNumber42 == 42:
            print('\n *** Invalid TIFF image. ***\n')
            return
        # Location of the first IFD
        IFDloc = struct.unpack('I', fid.read(4))[0]
        self.firstIFD = IFDloc
        # Go to first IFD
        fid.seek(IFDloc)
        # Get the number of IFD entries.
        NrEntries = struct.unpack('H', fid.read(2))[0]
        # Read entries
        header = {}
        for n in range(NrEntries):
            """
            Each IFD entry is 12 bytes long. An IFD entry is constructed as:
            (tag, data type, data count, value or pointer to value)
            """
            IFD = struct.unpack('HHII', fid.read(12))
            IFDtag = IFD[0]
            if IFDtag in self.IFDtags:
                TagDescription = self.IFDtags[IFDtag]
                if IFD[2] == 1:
                    #setattr(self, TagDescription, IFD[3])
                    header[TagDescription] = IFD[3]
                else:
                    #setattr(self, TagDescription, IFD)
                    header[TagDescription] = IFD
        
        fid.close()

        return header
    
    def read(self, frames=-1):
        w = self.header['ImageWidth']
        h = self.header['ImageHeight']

        if frames > 0:
            nrFrames = frames
        else:
            nrFrames = self.NrFrames

        NrPixels = w*h
        data = np.zeros((h, w, nrFrames))

        """ Image data locations"""
        #EOF = False
        #self.NrFrames = 0

        IFDloc = self.firstIFD

        fid = open(self.FilePath, 'rb')
        #while not EOF:
        for n in range(nrFrames):
            fid.seek(IFDloc)
            NrEntries = struct.unpack('H', fid.read(2))[0]

            IFDend = IFDloc + 2 + NrEntries*12
            dataLoc = IFDend + 4
            #self.DataLocs.append(IFDend + 4)

            fid.seek(dataLoc)
            data[:,:,n] = np.fromfile(fid, dtype=np.uint16, count=NrPixels).reshape(h, w)

            # Location of the next IFD
            fid.seek(IFDend)
            IFDloc = struct.unpack('I', fid.read(4))[0]
            
            #self.NrFrames += 1

            #if IFDloc == 0:
            #    EOF = True
        
        fid.close()
        #dt = float(self.log['Andor-ActualInterval-ms'])/1000.0
        #time = np.arange(self.NrFrames)*dt
        return Movie(data, header={})


    def _parse_log(self):
        fid = open(self.FilePath, 'rb')
        fid.seek(self.header['LogLoc'][3])
        LogLength = self.header['LogLoc'][2]
        Log = struct.unpack(str(LogLength)+'c', fid.read(LogLength))
        fid.close()
        
        LogStr = ''
        for s in Log:
            LogStr += s.decode(self.encoding)
        LogStr = LogStr.replace('null', 'None')
        LogStr = LogStr.replace('false', 'False')
        LogStr = LogStr.replace('true', 'True')
        
        #exec('LogDict = ' + LogStr)

        kvs = [kv.split(':') for kv in LogStr[1:-2].split(',')]

        logDict = {}
        for kv in kvs:
            logDict[kv[0][1:-1]] = kv[1][1:-1]
        
        return logDict



class IJMetaFile():
    """ Interface for ImageJ metadata as exported by uManager.

    Abstracts microscope frame properties as TimeSeries, e.g.
    time stamps, laser intensity etc.
    """

    def __init__(self, filePath):
        self.filePath = filePath
        self.summary, self.frameSeries = self.read()

    def read(self):
        fid = open(self.filePath)

        blocks = []
        for line in fid:
            if line[0] == '"':
                #block = line.split(':')[0].strip().strip('"')
                blocks.append({})
            elif line[0] == ' ':
                key = line.split(':')[0].strip().strip('"')
                arg = line.split(':')[1].strip()[:-1].strip('"')
                blocks[-1][key] = arg

        fid.close()

        summary = blocks[0]
        frames = blocks[1:]

        frameSeries = {}
        keys = frames[0].keys()
        for key in keys:
            frameSeries[key] = [m[key] for m in frames]

        return summary, frameSeries

    def get_timeseries(self, prop, dtype):
        
        if self.frameSeries.has_key(prop):
            meta = np.array(self.frameSeries[prop], dtype=dtype)
        else:
            meta = np.array([])

        return meta



class PmaFile():
    """ Interface for legacy binary files .pma
    
    Filetype for smFRET movies once created in TJ Ha lab and still used by some.
    """

    PixelBytes = 2
    InfoBytes = 4

    def __init__(self, FilePath):
        self.FilePath = FilePath

        self.header = self.get_header()
        
        #FrameSize = self.get_frame_size()
        #self.ImageWidth = FrameSize[1]
        #self.ImageHeight = FrameSize[0]

        #self.NrFrames = self.get_nr_frames()

        #self.load_frames()
        
    def get_frame_size(self):
        fid = open(self.FilePath, "rb")
        
        # The first 2x2 bytes represent the dimensions of the frame as uint16
        imsize_bin = fid.read(self.InfoBytes)
        imsize = struct.unpack('2H', imsize_bin)

        fid.close()
        return imsize

    def get_nr_frames(self, w, h):
        FileBytes = os.stat(self.FilePath)[stat.ST_SIZE]
        #FrameBytes = (self.ImageHeight*self.ImageWidth + 1)*self.PixelBytes
        FrameBytes = (h*w + 1)*self.PixelBytes
        return (FileBytes - self.InfoBytes)/FrameBytes

    def get_header(self):
        fileName = os.path.split(self.FilePath)[-1]

        FrameSize = self.get_frame_size()
        ImageWidth = FrameSize[1]
        ImageHeight = FrameSize[0]
        NrFrames = self.get_nr_frames(ImageWidth, ImageHeight)

        header = {
                'FILE': fileName,
                'WIDTH': ImageWidth,
                'HEIGHT': ImageHeight,
                'NR_FRAMES': NrFrames}
        return header

    def read(self, frames=-1):
        h = self.header['HEIGHT']
        w = self.header['WIDTH']
        
        NrFramePixels = h*w
        if frames > 0:
            NrFrames = frames
        else:
            NrFrames = self.header['NR_FRAMES']

        # Set the file pointer to the start of the first frame
        fid = open(self.FilePath, 'rb')
        FrameLoc = self.InfoBytes
        fid.seek(FrameLoc)

        # Read the frame data
        DataType = np.dtype([('frame', 'H'), ('data', 'H', NrFramePixels)])
        FrameArray = np.fromfile(fid, dtype=DataType, count=NrFrames)

        fid.close()
        
        self.FrameIndices = list(FrameArray['frame'])
        FrameArray = np.transpose(np.array(FrameArray['data']))
        FrameArray = np.reshape(FrameArray, (h, w, NrFrames))
        
        #fm.import_array(FrameArray)
        #fm.frames = FrameArray
        #dt = 0.1
        #time = np.arange(NrFrames)*dt
        #FrameArray *= 25.0

        header = self.get_header()
        #header['STRING'] = make_header_str(header)

        return Movie(FrameArray, header=header)


""" ===========================================================================
= TIMESERIES
===============================================================================

TimeSeries and associated metadata are stored in a custom plain text format,
with extension .trc for timeseries, .spc for metadata, and .sel for trace
selections. === These need to be replaced by hdf5 format ===

"""

class TraceFile():
    """ Interfacing with files storing TimeSeries data.
    
    Traces are stored in plain text files, in which each line, except for the first line, represents a 
    measurement from a single molecule. Lines can contain multiple blocks of data, in which each block 
    represents a different measurement or derived signal, e.g. donor signal, acceptor signal, calculated
    fret etc. The first line in the file contains descriptions of the blocks.
    """
    
    def __init__(self, filePath):
        if isinstance(filePath, TraceFile):
            self.filePath = filePath.filePath
        else:
            self.filePath = filePath

    def read_header(self):
        if self.exists():
            
            hstr = []
            f = open(self.filePath, 'r')
            for line in f:
                if line[0]!='#':
                    break
                line = line.rstrip('\n').strip()
                hstr.append(line)
            f.close()
        
            hstr = ' '.join(hstr)
            return self._parse_header_str(hstr)
        else:
            return {}
    
    @staticmethod
    def _parse_header_str(hstr):
        dct = {}
        for s in hstr.split('#')[1:]:
            key = s.split(':')[0].strip()
            val = s.split(':')[1].strip()

            keyFlag = key[-3:]
            val = val.strip().split()

            if len(val) == 1:
                val = val[0]

            dct[key] = val
        
        return dct

    @staticmethod
    def _make_header_str(dct):
        keys = list(dct.keys())
        values = list(dct.values())

        values = [' '.join(list(val)) for val in values]

        order = np.argsort(keys)
        hstr = ['# %s: %s'%(keys[c], values[c]) for c in order]

        return '\r\n'.join(hstr)

    def read(self):
        if os.path.isfile(self.filePath):
            blocks = self.read_header()['blocks']
            nrBlocks = len(blocks)
            
            hstr = []
            f = open(self.filePath, 'r')
            for line in f:
                if line[0]!='#':
                    line = line.rstrip('\n').strip().split()
                    line = np.array([float(l) for l in line])
                    hstr.append(line)
            f.close()

            dataSizes = [h.size for h in hstr]
            if len(hstr) == 0:
                return [[]], []

            maxSize = max(dataSizes)
            padding = (maxSize - np.array(dataSizes))/nrBlocks
            
            padding = np.array(padding, dtype=int)
            
            blockData = [np.hsplit(d, nrBlocks) for d in hstr]

            paddedData = [[np.hstack((d, np.nan*np.ones(padding[n]))) 
                for d in blockData[n]] for n in range(len(blockData))]
            
            data = []
            for n in range(nrBlocks):
                stackedData = np.vstack([d[n] for d in paddedData])
                data.append(stackedData)

            return [TimeSeries(d) for d in data], blocks

        else:
            return [], []

    def exists(self):
        return os.path.isfile(self.filePath)

    def clear(self):
        f = open(self.filePath, 'w')
        f.write('')
        f.close()

    def append(self, ts, descr):
        header = self.read_header()
        if self.exists():
            #if header.has_key('blocks'):
            if 'blocks' in header:
                #if cmp(descr, header['blocks'])==0:
                if len(descr)==len(header['blocks']):
                    self.write(ts, descr, append=True)
                else:
                    assert (len(descr)==0), 'File headers do not match.'
            else:
                if (len(header.keys())>0) & (len(descr)>0):
                    raise ValueError('Headers do not match while trying to append to file.')
                else:
                    self.write(ts, descr)
        else:
            self.write(ts, descr)

    def write(self, ts, descr=None, append=False):

        if isinstance(ts, list):
            
            # Check that all blocks have the same number of traces
            nrTraces = np.array([d.get_nr_traces() for d in ts])
            equalNrTraces = np.sum((nrTraces == nrTraces[0])) == len(ts)
            
            # Collect traces per molecule
            if equalNrTraces:
                if nrTraces[0] > 1:
                    data = np.vstack([np.hstack([d.data[n,:] for d in ts]) 
                        for n in range(nrTraces[0])])
                else:
                    data = np.hstack([d.data for d in ts])

            else:
                raise ValueError('Blocks need equal number of traces when writing to file.')

            if isinstance(descr, list):

                # Check that number of descriptions matches number of blocks
                equalNrDescr = len(descr) == len(ts)

                if not equalNrDescr:
                    raise ValueError('Each blocks needs a description when writing to file.')
            else:
                raise ValueError('Each blocks needs a description when writing to file.')

            nrTraces = nrTraces[0]

        else:
            data = ts.data
            descr = list(descr)[0]
            nrTraces = ts.get_nr_traces()
        
        #cdate = datetime.datetime.now().strftime('%Y-%m-%d')
        #header = {'blocks': descr, 'timestamp':[cdate]}
        header = {'blocks': descr}
        headerStr = self._make_header_str(header)

        #np.savetxt(self.filePath, data, header=headerStr, fmt='%.3f', *args, **kwargs)
        
        if append:
            mode = 'a'
        else:
            mode = 'w'
        f = open(self.filePath, mode)

        if mode == 'w':
            f.write(headerStr)
            f.write('\n')
        
        if nrTraces > 1:
            for n in range(nrTraces):
                [f.write('%.3f '%d) for d in data[n,:]]
                f.write('\n')
        elif nrTraces == 1:
            data = np.reshape(data, data.size)
            [f.write('%.3f '%d) for d in data]
            f.write('\n')

        f.close()
        
        #if not append:
        #    print 'Data saved to %s'%os.path.split(self.filePath)[-1]

    def write_block(self, ts, descr):
        data, blocks = self.read()
        if descr in blocks:
            data[blocks.index(descr)] = ts
        else:
            data.append(ts)
            blocks.append(descr)
        self.write(data, blocks)




class SpecsFile(object):
    """ Interfacing with files storing metadata. 
    
    """

    def __init__(self, filePath):
        self.filePath = filePath

    def read(self):
        if os.path.isfile(self.filePath):
            try:
                f = open(self.filePath, 'r')
                headerStr = f.readline()
                f.close()

                header = headerStr.lstrip('#').strip().split()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    data = np.loadtxt(self.filePath, skiprows=1)
                specs = SmfSpecs(data, header)
            except:
                specs = SmfSpecs()
        else:
            specs = SmfSpecs()
            
        return specs

    def exists(self):
        return os.path.isfile(self.filePath)
    
    def clear(self):
        f = open(self.filePath, 'w')
        f.write('')
        f.close()

    def write(self, specs):
        headerStr = ' '.join(specs.descr)
        data = specs.data
        if len(data.shape) < 2:
            data = np.reshape(data, [1,data.size])

        if isinstance(specs, SmfSelect):
            fmt = '%d'
        else:
            fmt = '%.3f'
        
        return np.savetxt(self.filePath, data, fmt=fmt, header=headerStr)
    
    def append(self, specs):
        if self.exists():
            content = self.read()
            content.append(specs)
        else:
            content = specs

        self.write(content)



class SelectFile(SpecsFile):
    """ Interfacing with files storing trace selections. 
    
    """

    def read():
        pass

    def write():
        pass




""" ===========================================================================
= DATASETS
===============================================================================
"""


class AnalysisFile():
    """ Interfacing with files storing datasets.

    Datasets are stored in plain text files, and specified as relative paths to
    data folders. They can be grouped under a common header describing eg which
    sample the data came from.
    """

    def __init__(self, filePath):
        self.filePath = filePath

    def read(self):
        f = open(self.filePath, 'r')
        groups = []
        data = []
        header = []
        for line in f:
            if line[0]=='#':
                header.append(line.lstrip('#').strip())
            elif line[0]=='<':
                # New group
                grpName = line.strip()[1:-1]
                groups.append(grpName)
                data.append([[]])
            else:
                content = line.strip('\n')
                if content=='':
                    data[-1].append([])
                else:
                    data[-1][-1].append(content)
        
        f.close()

        header = AnalysisFile._parse_header_lines(header)

        return groups, data, header

    def write(self, groups, data, header={}):
        f = open(self.filePath, 'w')

        header = AnalysisFile._build_header_lines(header)

        [f.write('# %s\n'%h) for h in header]

        for n in range(len(groups)):
            key = groups[n]
            grpData = data[n]

            f.write('<'+key+'>\n')

            #grpData = data[key]
            dataStrings = []
            for d in grpData:
                if not isinstance(d, np.ndarray):
                    d = np.array([d])
                dataStrings.append(('%.3f '*d.size)%tuple(d)+'\n')
            f.writelines(dataStrings)
        
        f.close()

    
    @staticmethod
    def _parse_header_lines(headerLines):
        headerDict = {}
        for h in headerLines:
            info = h.split(':')
            headerDict[info[0].strip()] = info[1].strip()
        return headerDict

    @staticmethod
    def _build_header_lines(headerDict):
        headerLines = []
        for key in headerDict:
            headerLines.append('%s: %s'%(key, headerDict[key]))
        return headerLines

