""" Movie analysis routines.
"""


# System packages
import os, glob
import numpy as np
from scipy.ndimage import generic_filter, filters
from scipy.ndimage.morphology import binary_erosion, binary_dilation, \
        grey_erosion
from scipy.ndimage.measurements import label, find_objects, \
        center_of_mass, maximum_position
from scipy.ndimage.interpolation import rotate, geometric_transform
from scipy.ndimage import median_filter

# Program modules
from smfproc.source.io.datatypes import TimeSeries, ImageCoordinates, Frame
from smfproc.source.lib.math import nanmedian
from smfproc.source.lib.tools import list_apply, list_diff_apply, list_logical_or, \
        list_logical_select


def scale_fluor_intensity(don, acc, crds, scaleMatrix):
    #row, col = crds.get_row_col()
    row = crds[:,0]
    col = crds[:,1]
    nrTraces = don.get_nr_traces()
    for n in range(nrTraces):
        scale = scaleMatrix[row[n], col[n]]
        #scale[scale==0] = 1.0
        don.data[n,:] = don.data[n,:]/scale
        acc.data[n,:] = acc.data[n,:]/scale


def get_image_background(im, k=30):
    #return median_filter(im, kernel, mode='mirror')
    im = np.reshape(im, im.shape[:2])
    nrRows, nrCols = im.shape
    nrHorBlocks = int(nrCols/k)+1
    nrVertBlocks = int(nrRows/k)+1
    imMed = np.zeros(im.shape)
    for c in range(nrHorBlocks):
        for r in range(nrVertBlocks):
            values = im[r*k:(r+1)*k, c*k:(c+1)*k]
            imMed[r*k:(r+1)*k, c*k:(c+1)*k] = np.nanmean(values) #nanmedian(values)
    return imMed



def analyze_movie(movies, warpMatrix, movieSettings):
    #print('%s - Detecting molecules...' % mFile.filePath)
    sett = movieSettings

    #trcs, crds, rdispl, cdispl = get_traces_from_movie(movies, 
    donSignals, accSignals, donCrds, accCrds, donCOM, accCOM = \
            get_traces_from_movie(movies,
                    warpMatrix,
                    sett['SIGNAL_FRAMES'],
                    sett['PIXEL_SNR'],
                    sett['MOLECULE_RADIUS'],
                    sett['BACKGROUND_RADIUS'])
    
    trcs = [
            TimeSeries(donSignals),
            TimeSeries(accSignals),
            TimeSeries(donCOM),
            TimeSeries(accCOM)]

    crds = [donCrds, accCrds]
    
    return trcs, crds


def get_traces_from_movie(movs, 
        warpMatrix,
        nrSignalFrames,
        snr,
        signalRadius,
        bckGrndRadius,
        detectLeft=True,
        detectRight=True):
    
    # Detect signals
    donMov = movs[0].frames
    accMov = movs[1].frames
    
    donMaxIm = get_max_fluorescence(donMov, nrSignalFrames)
    accMaxIm = get_max_fluorescence(accMov, nrSignalFrames)

    donCrds, donMask = detect_molecules(donMaxIm, snr)
    accCrds, accMask = detect_molecules(accMaxIm, snr)
    
    # Register right and left side molecules
    donCrds, accCrds = register_crds(donCrds, accCrds, signalRadius, 
            warpMatrix, directMapOnly=True)
    
    # Remove border molecules
    sr = signalRadius
    borders = np.array((sr, donMov.shape[0]-sr, sr, donMov.shape[1]-sr))
    donCrds, accCrds = remove_border_crds(donCrds, accCrds, borders)

    # Create signal masks
    signalROI = get_circular_roi(signalRadius)
    bckGrndROI = get_circular_roi(bckGrndRadius)
    crds2im(donMask, donCrds, signalROI, True)
    crds2im(accMask, accCrds, signalROI, True) 

    # Get traces
    donSignals, donCOMs = get_signals(
            donMov, donMask, donCrds, signalROI, bckGrndROI)
    accSignals, accCOMs = get_signals(
            accMov, accMask, accCrds, signalROI, bckGrndROI)
 
    return donSignals, accSignals, donCrds, accCrds, donCOMs, accCOMs
    

def get_signals(mov, mask, crds, sroi, broi):
    nrCrds = crds.shape[0]
    nrFrames = mov.shape[2]

    signals = np.zeros((nrCrds, nrFrames))
    coms = np.zeros((nrCrds, nrFrames))

    for n in range(nrCrds):
        sroic = clip_roi(sroi, crds[n,:], mask.shape)
        broic = clip_roi(broi, crds[n,:], mask.shape)

        signal, com = get_signal(mov, mask, sroic, broic) 
        signals[n,:] = signal
        coms[n,:] = com

    return signals, coms


def get_signal(mov, mask, sroi, broi):
    
    signal = get_trace(mov, np.ones(mask.shape)==True, sroi) 
    bckGrnd = get_trace(mov, mask==False, broi)
    com = get_com_displacement(sroi, signal)
    
    signal = np.sum(signal, axis=0) 
    #roiCrds = get_crd_roi(crd, sroi)

    ROIarea = sroi.shape[0]
    bckGrndPerFrame = np.median(bckGrnd, axis=0)
    #bckGrndPerFrame = np.reshape(bckGrndPerFrame, [1, len(bckGrndPerFrame)])
    #signal -= np.repeat(bckGrndPerFrame, ROIarea, 0)
    #trace = np.sum(signal, axis=0)

    

    trace = signal - ROIarea*bckGrndPerFrame
    
    #plt.figure()
    #plt.plot(trace)
    #plt.show()
    #exit()

    return trace, com


def clip_roi(roi, crd, shp):
    r = crd[0]+roi[:,0]
    c = crd[1]+roi[:,1]

    crit = ((r >= 0) & (r < shp[0])) & ((c >= 0) & (c < shp[1]))
    roi = np.vstack((np.array(r[crit]), np.array(c[crit])))
    
    return np.transpose(roi)


def get_com_displacement(crds, traces):
    traces = np.array(traces, dtype=float)
    
    nrTraces = traces.shape[0]
    length = traces.shape[1]

    rows = np.reshape(crds[:,0], [nrTraces, 1])
    cols = np.reshape(crds[:,1], [nrTraces, 1])

    rowCrds = np.repeat(rows, length, 1)
    colCrds = np.repeat(cols, length, 1)

    center = np.mean(crds, axis=0)
    rowCrds = rowCrds - center[0]
    colCrds = colCrds - center[1]

    tracesSum = np.sum(traces, axis=0)
    rcom = np.sum(rowCrds*traces, axis=0)/tracesSum
    ccom = np.sum(colCrds*traces, axis=0)/tracesSum

    displ = np.sqrt(rcom**2 + ccom**2)

    return displ


def get_trace(mov, mask, roi):
    #r = crd[0]+roi[:,0]
    #c = crd[1]+roi[:,1]
    r = roi[:,0]
    c = roi[:,1]
    trace = mov[r, c, :]
    crit = mask[r, c]
    return trace[crit,:]


def crds2im(im, crds, roi, value):
    nrCrds = crds.shape[0]
    for n in range(nrCrds):
        nroi = get_crd_roi(crds[n,:], roi)
        im[nroi[:,0], nroi[:,1]] = value


def get_crd_roi(crd, roi):
    nroi = np.array(roi)
    nroi[:,0] += crd[0]
    nroi[:,1] += crd[1]
    return nroi


def get_circular_roi(radius):
    crds = []
    for r in range(-radius, radius+1):
        for c in range(-radius, radius+1):
            d = np.sqrt(c**2 + r**2)
            if d <= radius:
                crds.append([r, c])
    return np.array(crds)


def register_crds(donCrds, accCrds, overlap, warpMatrix, directMapOnly=False):
    accCrdsInDonFrame = shift_crds(accCrds, warpMatrix)

    D = get_crd_distances(donCrds, accCrdsInDonFrame)
    
    crit = D < overlap+1
    hitDon, hitAcc = np.where(crit)

    donMapped = donCrds[hitDon, :]
    accMapped = accCrds[hitAcc, :]

    accNotMapped = accCrds[np.sum(crit, axis=0)==0, :]
    donNotMapped = donCrds[np.sum(crit, axis=1)==0, :]

    accCrdsInDonFrame = shift_crds(accNotMapped, warpMatrix)
    donCrdsInAccFrame = shift_crds(donNotMapped, warpMatrix, -1)

    if directMapOnly:
        donCrds = donMapped
        accCrds = accMapped
    else:
        donCrds = np.vstack((donMapped, donNotMapped, accCrdsInDonFrame))
        accCrds = np.vstack((accMapped, donCrdsInAccFrame, accNotMapped))
 
    return donCrds, accCrds


def get_crd_distances(don, acc):
    # Make an nxm matrix of don-acc distance
    nrDonCrds = don.shape[0]
    nrAccCrds = acc.shape[0]

    D = np.zeros((nrDonCrds, nrAccCrds))
    
    for n in range(nrDonCrds):
        for m in range(nrAccCrds):
            donRow = don[n,0]
            donCol = don[n,1]
            accRow = acc[m,0]
            accCol = acc[m,1]
            dist = np.sqrt((donRow-accRow)**2 + (donCol-accCol)**2)
            D[n,m] = dist

    return D


def remove_border_crds(donCrds, accCrds, borders):
    donCrit = is_within_border(donCrds, borders)==False
    accCrit = is_within_border(accCrds, borders)==False
    crit = np.logical_and(donCrit, accCrit)
    return donCrds[crit,:], accCrds[crit,:]


def is_within_border(crd, b):
    withinBorder = list_logical_or([
        crd[:,0]<b[0], 
        crd[:,0]>=b[1], 
        crd[:,1]<b[2], 
        crd[:,1]>=b[3]])
    return withinBorder[-1:][0]

    
def shift_crds(crds, M, norm=1.0):
    outCrd = np.array(crds)
    for i in range(crds.shape[0]):
        dx = M[0][int(crds[i][0]), int(crds[i][1])]
        dy = M[1][int(crds[i][0]), int(crds[i][1])]
        outCrd[i,0] = crds[i][0]+norm*dy
        outCrd[i,1] = crds[i][1]+norm*dx
    return outCrd

    
def get_max_fluorescence(data, res=5):
    """ Split the movie data into blocks of 'res' frames, then average those.
    Project the maximal signal for each pixel in the blocked sequence.

    input arguments:
    data -> movie frames [numpy array]
    res -> block size [int]

    output:
    projected signal [numpy array]
    """

    nrFrames = int(data.shape[2]/res)
    blockedIm = []
    for m in range(nrFrames):
        blockedIm.append(np.mean(data[:,:,res*m:res*(m+1)], axis=2))

    A = np.array(blockedIm)
    B = np.max(A, axis=0)

    return B


def detect_molecules(image, snr):

    # Remove dead/unilluminated pixels
    globalMin = np.min(image)
    globalMed = np.median(image)
    deadCriterium = image < (globalMin + 0.3*(globalMed-globalMin))
    image[deadCriterium] = np.nan

    # Make signal mask
    localMed, localStd = get_image_bgr(image, 50)
    globalStd = np.median(localStd)
    
    image -= localMed
    image[deadCriterium] = 0
    
    mask = image > snr*globalStd
    #mask = binary_dilation(binary_erosion(mask, iterations=1), iterations=1)

    # Get molecule coordinates
    signal = image*mask
    crds = get_local_maxima(signal, 5)
    
    return np.array(crds, dtype=int), mask


def get_image_bgr(im, res): 
    im = np.squeeze(im)

    nrRows, nrCols = im.shape

    k = res

    nrHorBlocks = int(nrCols/k)+1
    nrVertBlocks = int(nrRows/k)+1

    colEdges = np.array(np.round(np.linspace(0, nrCols, nrHorBlocks)), dtype=int)
    rowEdges = np.array(np.round(np.linspace(0, nrRows, nrVertBlocks)), dtype=int)

    imMed = np.zeros(im.shape)
    imStd = np.zeros(im.shape)
    for c in range(nrHorBlocks-1):
        for r in range(nrVertBlocks-1):
            values = im[rowEdges[r]:rowEdges[r+1], colEdges[c]:colEdges[c+1]]

            imMed[rowEdges[r]:rowEdges[r+1], colEdges[c]:colEdges[c+1]] = \
                    np.nanmean(values)
            imStd[rowEdges[r]:rowEdges[r+1], colEdges[c]:colEdges[c+1]] = \
                    np.nanstd(values)
    
    return imMed, imStd


def get_local_maxima(image, neighborhood):
    
    data_max = filters.maximum_filter(image, neighborhood)
    maxima = (image == data_max)
    # maxima now only has single pixels at the peaks
    # if two maxima are too close then this procedure picks the brightest one - is that what
    # we want?

    labeled, numObjects = label(maxima)
    # labeled objects are only one pixel size

    # NB: This gives some warning
    # NB: this doesn't seem right. taking the com of one-pixel objects. 
    # use find_objects instead?
    xy = np.array(center_of_mass(image, labeled, range(1, numObjects+1)))
    print('(Warning in get_local_maxima in movie analysis)')
    xy = xy[np.isnan(xy)==False]
    xy = np.reshape(xy, [xy.size/2, 2])

    return xy


