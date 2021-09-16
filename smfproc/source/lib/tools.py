""" Some general manipulations for lists
"""

import numpy as np

list_apply = lambda L, f, *a, **kwa: [f(l, *a, **kwa) for l in L]
list_diff_apply = lambda f, *A: \
        [f(*[a[i] for a in A]) for i in range(len(A[0]))]

list_logical_select = lambda L, c: [L[i] for i in np.nonzero(c)[0]]

lor = np.logical_or
def list_logical_or(L):
    for i in range(1,len(L)):
        L[i] = np.logical_or(L[i], L[i-1])
    return L

