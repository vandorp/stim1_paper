import os, glob
import numpy as np

from scipy.ndimage import generic_filter, filters
from scipy.ndimage.morphology import binary_erosion, binary_dilation, \
        grey_erosion
from scipy.ndimage.measurements import label, find_objects, \
        center_of_mass, maximum_position
from scipy.ndimage.interpolation import rotate, geometric_transform
from scipy.ndimage import median_filter

from ..core.datatypes import TimeSeries
from ..core.analysis import nanmedian

from misc import list_apply, list_diff_apply, list_logical_or, list_logical_select
from datatypes import ImageCoordinates, Frame



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


""" --- MOVIE ANALYSIS ----------------------------------------------------- """


def analyze_movie(movies, warpMatrix, movieSettings):
    #print('%s - Detecting molecules...' % mFile.filePath)
    sett = movieSettings

    trcs, crds, rdispl, cdispl = get_traces_from_movie(movies, 
            warpMatrix,
            sett['SIGNAL_FRAMES'],
            sett['PIXEL_SNR'],
            sett['MOLECULE_RADIUS'],
            sett['BACKGROUND_RADIUS'])
    
    displacement = [rdispl[0], cdispl[0]]
    
    if not sett['DONOR_LEFT']:
        trcs = trcs[::-1]
        crds = crds[::-1]
        displacement = displacement[::-1]

    trcs.extend(displacement)
    
    return trcs, crds


def get_traces_from_movie(movs, warpMatrix,
        nrSignalFrames=5,
        snr=1,
        signalRadius=4,
        bckgrndRadius=12,
        detectLeft=True,
        detectRight=True):

    # Get molecule coordinates from both sides
    signalFrames = range(nrSignalFrames)
    allCrds = list_apply(movs, detect_molecules, signalFrames, snr)

    # Map right side molecules to left side molecules
    crds = map_sm_crds_warp_matrix(
            allCrds, signalRadius, warpMatrix,
            detectleft=detectLeft,
            detectright=detectRight)

    # Get trace time series
    trcs = list_diff_apply(get_traces, movs, crds, 
            [signalRadius]*2, [bckgrndRadius]*2, allCrds)

    signal = [t[0] for t in trcs]
    rdispl = [t[1] for t in trcs]
    cdispl = [t[2] for t in trcs]

    return signal, crds, rdispl, cdispl


def map_sm_crds_warp_matrix(crds, overlap, warpMatrix, 
        detectleft=True, detectright=True,
        border=(0,0,0,0)):

    border_width = np.array(border) + overlap

    # shift acceptor crds to donor frame
    shift_crds_matrix(crds[1], warpMatrix)

    # merge donor and acceptor coordinates
    if (detectleft & detectright):
        crd = _merge_crds(crds)
    elif (detectright & (not detectleft)):
        crd = crds[1]
    elif (detectleft & (not detectright)):
        crd = crds[0]
    else:
        crd = crds[0]
        crd.crds = []
    
    # remove border objects
    _remove_border_crds(crd, border_width)
    shift_crds_matrix(crd, [-warpMatrix[0], -warpMatrix[1]])
    _remove_border_crds(crd, border_width)
    
    # detect molecules
    im = _crds2image(crd)
    accCrds = _get_object_crds(im, radius=1)

    # remove overlapping molecules
    accCrds = remove_overlapping_crds(accCrds, overlap)

    donCrds = ImageCoordinates(accCrds)
    shift_crds_matrix(donCrds, warpMatrix)

    return [donCrds, accCrds]

def remove_overlapping_crds(crds, radius):
    D = crds.get_distances()
    removeCrit = D < (2*radius+1)
    crit = np.sum(removeCrit, 0) > 1
    return crds.select(crit==False)

def shift_crds_matrix(crd, M):
    crdl = crd.get()
    for i in range(len(crdl)):
        dx = M[0][int(crdl[i][0]), int(crdl[i][1])]
        dy = M[1][int(crdl[i][0]), int(crdl[i][1])]
        crdl[i] = (crdl[i][0]+dy, crdl[i][1]+dx)

"""
def get_image_xtalk(don, acc):
    don = np.reshape(don.data, don.data.size)
    acc = np.reshape(acc.data, acc.data.size)
    slope = acc/don
    
    xtalk = []
    for nrBins in np.arange(20,60,5):
        cnt, binEdges = np.histogram(slope, bins=nrBins, range=(0,1))
        bins = (binEdges[0:-1]+binEdges[1:])/2.0
        xtalk.append(bins[np.argmax(cnt)])
    
    return np.median(xtalk)
"""
def get_traces(mov, crds, mol_rad, bgr_rad, bckGrndCrds):
    
    signalMask = _crds2image(crds, radius=mol_rad)
    bckGrndSignalMask = _crds2image(bckGrndCrds, radius=mol_rad)

    meanFrame = np.mean(mov.frames, axis=2)
    imageMask = meanFrame > np.std(meanFrame)
    
    nrMolecules = len(crds.get())
    nrFrames = mov.nr_frames()
    
    trcs = np.zeros((nrMolecules, nrFrames))
    rcom = np.zeros((nrMolecules, nrFrames))
    ccom = np.zeros((nrMolecules, nrFrames))
    for n in range(nrMolecules):

        # Get signal
        signalROI = get_crd_roi(crds.get_crd(n), mol_rad)
        src = signalROI.get_row_col()
        signal = mov.get()[src[0], src[1], :]
        signalSum = np.sum(signal, axis=0)
        signalArea = len(signalROI.get())
        
        # Get background
        backgrndMaskROI = _crds2image(crds.get_crd(n), bgr_rad)
        backgrndMask = imageMask*(bckGrndSignalMask.get()==0)*backgrndMaskROI.get()
        brc = np.where(backgrndMask)
        background = np.median(mov.get()[brc[0], brc[1], :], axis=0)
        
        trcs[n,:] = signalSum - background*signalArea
 
        bckgrnd = np.transpose(np.reshape(np.repeat(background*signalArea, signalArea), 
            [nrFrames, signalArea]))

        #trcSignal = signal-bckgrnd
        trcSignal = signal
        trcSignalSum = np.sum(trcSignal, axis=0)
        
        # Get signal center of mass
        centerCrds = crds.get_crd(n).crds[0]

        # center around molecule coordinate
        rrc = np.reshape(np.repeat(src[0]-centerCrds[0], nrFrames), 
                [signalArea, nrFrames])
        crc = np.reshape(np.repeat(src[1]-centerCrds[1], nrFrames), 
                [signalArea, nrFrames])
        
        rcom[n,:] = np.array(np.sum(rrc*trcSignal, axis=0)/trcSignalSum, dtype=float)
        ccom[n,:] = np.array(np.sum(crc*trcSignal, axis=0)/trcSignalSum, dtype=float)

        # This is a hack to get kind of pixel units
        rcom[n,:] *= np.sqrt(signalArea)
        ccom[n,:] *= np.sqrt(signalArea)

    #return TimeSeries(trcs[:,:500], mov.time[:500]), TimeSeries(rcom[:,:500], mov.time[:500]), TimeSeries(ccom[:,:500], mov.time[:500])
    
    return TimeSeries(trcs), TimeSeries(rcom), TimeSeries(ccom)


import matplotlib.pyplot as plt
def detect_molecules(mov, frames, snr):

    # Get average images for analysis
    signalFrame = mov.get_average_frame(frames)
    result = generic_filter(signalFrame.get(), nanmedian, size=50)
    
    # make binary images
    mask_signals(signalFrame, result, snr)
    #create_mask_for_object_tagging(signalFrame, result, snr)

    #xy = get_local_maxima(signalFrame, 5)
    
    #plt.figure()
    #plt.imshow(signalFrame.get())
    #plt.plot(xy[:, 1], xy[:, 0], 'ro')
    #plt.show()
    #exit()

    return _get_object_crds_max(signalFrame)



def get_local_maxima(im, neighborhood):
    data = im.get()
    data_max = filters.maximum_filter(data, neighborhood)
    maxima = (data == data_max)

    labeled, num_objects = label(maxima)

    plt.figure()
    plt.imshow(labeled)
    xy = np.array(center_of_mass(data, labeled, range(1, num_objects+1)))

    return xy


def create_mask_for_object_tagging(im, bckGrnd, snr=0, imThresh=0.25, sigMaskThresh=0.2):
    pix = im.get()

    offset = np.min(pix)

    pix -= offset
    bckGrnd -= offset

    crit = pix < imThresh*np.median(pix)

    pix2 = np.array(pix)
    pix2[crit] = np.nan

    result = bckGrnd

    mask = pix > (result + snr*np.nanstd(pix2))
    mask = binary_dilation(binary_erosion(mask, iterations=1), iterations=1)
    
    signal = pix*mask

    im.set(signal)
    
    
def mask_signals(im, bckGrnd, snr=0, imThresh=0.25, sigMaskThresh=0.2):
    pix = im.get()

    offset = np.min(pix)

    pix -= offset
    bckGrnd -= offset

    crit = pix < imThresh*np.median(pix)

    pix2 = np.array(pix)
    pix2[crit] = np.nan

    result = bckGrnd

    mask = pix > (result + snr*np.nanstd(pix2))
    mask = binary_dilation(binary_erosion(mask, iterations=1))
    
    signal = pix*mask
    
    # Refine mask
    objs = find_objects(label(mask)[0])
    for obj in objs:
        isolatedSignal = signal[obj]

        subMask = isolatedSignal > sigMaskThresh*np.max(isolatedSignal)
        nrElements = np.sum(subMask)
        while nrElements > 20:
            subMask = binary_erosion(subMask)
            nrElements = np.sum(subMask)

        mask[obj] = subMask
    
    im.set(mask*pix)

def _get_object_crds_max(im):
    pix = im.get()
    lab = label(pix>0)[0]
    nrObjs = np.max(lab)
    center = maximum_position(pix, labels=lab, index=range(1,nrObjs+1))
    return ImageCoordinates(center, im.shape())

def _get_object_crds(im, radius=0):
    pix = im.get()
    if radius > 0:
        pix = binary_dilation(pix, iterations=radius)
    lab = label(pix)[0]
    nrObjs = np.max(lab)
    com = center_of_mass(pix, labels=lab, index=range(1,nrObjs+1))
    return ImageCoordinates(com, im.shape())


def shift_crds_matrix(crd, M):
    crdl = crd.get()
    for i in range(len(crdl)):
        dx = M[0][int(crdl[i][0]), int(crdl[i][1])]
        dy = M[1][int(crdl[i][0]), int(crdl[i][1])]
        crdl[i] = (crdl[i][0]+dy, crdl[i][1]+dx)

def _merge_crds(crds):
    mergedCrds = crds[0].get()+crds[1].get()
    return ImageCoordinates(mergedCrds, crds[0].imageShape)

def _remove_border_crds(crd, border):
    crdl = crd.get_row_col()
    b = np.array(border)+1
    w = crd.imageShape[1]
    h = crd.imageShape[0]
    within_border = list_logical_or( \
            [crdl[0]<b[2], crdl[0]>=(h-b[3]), crdl[1]<b[0], crdl[1]>=(w-b[1])])[-1:][0]

    selCrds = list_logical_select(crd.get(), within_border==False)
    crd.set(selCrds)

def get_crd_roi(center, radius=0):
    cnt = center.get()[0]
    shp = center.imageShape

    crds = []
    for r in range(-radius, radius+1):
        for c in range(-radius, radius+1):
            d = np.sqrt(c**2 + r**2)
            if d <= radius:
                row = r + cnt[0]
                col = c + cnt[1]
                if (row >= 0) & (row < shp[0]) & (col >= 0) & (col < shp[1]):
                    crds.append((row, col))
    return ImageCoordinates(crds, shp)

def _crds2image(crds, radius=0):
    im = np.zeros(crds.imageShape)
    for i in range(len(crds.get())):
        crd = crds.get_crd(i)
        im[crd.get()[0][0], crd.get()[0][1]] = 1
        if radius > 0:
            roi = get_crd_roi(crd, radius)
            rc = roi.get_row_col()
            im[rc[0], rc[1]] = 1
    return Frame(im)

