#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt('cond_aniani.out')
nrm = np.loadtxt('norm_aniani.out')

ndx = np.invert(np.isfinite(dat))
dat[ndx] = 0.
dat *= nrm

#ll = dat.copy()
ll = np.sum(dat, axis=1)

dat = np.loadtxt('cond_anicat.out')
nrm = np.loadtxt('norm_anicat.out')

ndx = np.invert(np.isfinite(dat))
dat[ndx] = 0.
dat *= nrm

#ll += dat.copy()
ll += np.sum(dat, axis=1)

dat = np.loadtxt('cond_catcat.out')
nrm = np.loadtxt('norm_catcat.out')

ndx = np.invert(np.isfinite(dat))
dat[ndx] = 0.
dat *= nrm

#ll += dat.copy()
ll += np.sum(dat, axis=1)

dat = np.loadtxt('cond_neinst.out')
nrm = np.loadtxt('norm_neinst.out')

ndx = np.invert(np.isfinite(dat))
dat[ndx] = 0.
dat *= nrm

ll += dat

#plt.plot(ll)
plt.plot(ll / nrm, 'x')

dat = np.loadtxt('cond_all.out')

plt.plot(dat)

plt.show()
