#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

vx = np.loadtxt('vx.out')
vy = np.loadtxt('vy.out')
vz = np.loadtxt('vz.out')

def correlate(a, b):
    fta = np.fft.rfft(a, n=2*len(a))
    ftb = np.fft.rfft(b, n=2*len(b))

    prd = fta*np.conj(ftb)

    mycor = np.fft.irfft(prd, n=2*len(prd))

    return mycor[:len(mycor)/2]

corx = None

for i in xrange(np.shape(vx)[1]):
    for j in xrange(i, np.shape(vx)[1]):

        forc = correlate(vx[:,i], vx[:,j])

        if i != j:
            forc *= 2

        #if i % 100 == 0:
        #    plt.plot(forc, label='%i-%i' %(i,j))

        if corx is None:
            corx = forc.copy()
        else:
            corx += forc

cory = None

for i in xrange(np.shape(vy)[1]):
    for j in xrange(i, np.shape(vy)[1]):

        forc = correlate(vy[:,i], vy[:,j])

        if i != j:
            forc *= 2

        #if i % 100 == 0:
        #    plt.plot(forc, label='%i-%i' %(i,j))

        if cory is None:
            cory = forc.copy()
        else:
            cory += forc

corz = None

for i in xrange(np.shape(vz)[1]):
    for j in xrange(i, np.shape(vz)[1]):

        forc = correlate(vz[:,i], vz[:,j])
        #if i % 100 == 0:
        #    plt.plot(forc, label='%i-%i' %(i,j))

        if i != j:
            forc *= 2

        if corz is None:
            corz = forc.copy()
        else:
            corz += forc

corr = (corx + cory + corz ) / np.shape(vx[1])

plt.plot(corr, linewidth=2.)

plt.legend()
plt.show()
