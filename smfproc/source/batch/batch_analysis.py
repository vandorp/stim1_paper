
# Python packages
import numpy as np
import matplotlib.pyplot as plt

# Program modules
from smfproc.source.io.filetypes import MovieFile


def load_movie_data(filePath, nrFrames=-1):
    movieFile = MovieFile(filePath)
    movie = movieFile.read(nrFrames)
    movs = movie.split()
    header = movieFile.header
    if nrFrames > 1:
        for mov in movs:
            mov.frames = np.mean(mov.frames, axis=2)
            mov.frames = np.reshape(mov.frames, 
                    [mov.frames.shape[0], mov.frames.shape[1],1])
    return movs, header


def get_background_illum(moviePath, resolution):
    movs,_ = load_movie_data(moviePath, 1)

    bckgr = get_image_background(movs[0].frames, k=resolution)
    
    bckgr = np.reshape(bckgr, bckgr.shape[:2])
    #bckgr = bckgr/np.max(bckgr)
    #bckgr = bckgr - 0.45
    #bckgr[bckgr<0] = 0

    return bckgr


def generate_laser_profile(moviePaths, resolution):
    """ Generate an average laser profile from a batch of movies. """

    bckgr = []
    for mPath in moviePaths:
        illum = get_background_illum(mPath, resolution)
        bckgr.append(illum)
    bckgr = np.array(bckgr)
    bckgr = np.mean(bckgr, axis=0)
    
    return bckgr
