from hmm import get_hmm, get_best_hmm

from smfproc.core.filetypes import TraceFile

import numpy as np

import matplotlib.pyplot as plt

data = TraceFile('selected.trc').read()

index = 35

don = data[0][0].get(index).data
acc = data[0][1].get(index).data
data = np.transpose(np.vstack((acc, don)))

#nStates = 3
#params = get_hmm(data, nStates)
params = get_best_hmm(data)

P = params.postFit
#N = nStates
N = P.shape[0]

#if N==1:
#    P = np.reshape(P, [1,P.size])

# make model
stateIndex = np.argmax(P, 0)
S = np.zeros((N, don.size))

for n in range(don.size):
    S[stateIndex[n], n] = 1

M = np.zeros(don.size)
for n in range(N):
    crit = S[n,:]==1
    M[crit] = np.mean(don[crit])

ax = plt.subplot(211)
plt.plot(don, 'g')
plt.plot(acc, 'r')
plt.plot(M, 'k')
ax.set_xlim([0, 510])

plt.subplot(212)
plt.plot(np.transpose(P))
plt.axis([0,510,-0.2,1.2])

plt.show()
