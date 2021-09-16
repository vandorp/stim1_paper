""" General calculations. """

import numpy as np

def nanmean(data):
    data = np.array(data)
    crit = np.isnan(data) == False
    data = data[crit]

    if data.size == 0:
        return np.nan
    else:
        return np.mean(data)

def nanmedian(d):
    m = np.array(d)
    m = m[np.isnan(m)==False]
    return np.median(m)

