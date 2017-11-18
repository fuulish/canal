#!/usr/bin/env python

import numpy as np

from numpy.random import randint
from numpy import arange
from numpy import vstack

nmol = 10
ntimesteps = 1000

a = arange(0, ntimesteps)

for i, c in zip(range(3), ['x', 'y', 'z']):
    #10 molecules with starting values between [0, 100)
    mols = randint(0, ntimesteps, (nmol))

    #just the progression of the coordinate
    xyzp = vstack([a]*nmol).T

    xyzp += mols

    np.savetxt('%scom.dat' %c, xyzp, fmt='%7.1f')

chg = open('charges.dat', 'w')

for n in range(nmol/2):
    chg.write(' 1.\n')

for n in range(nmol/2):
    chg.write('-1.\n')

chg.close()

cell = np.ones(ntimesteps)*ntimesteps*10
np.savetxt('cell.dat', cell)

print xyzp
