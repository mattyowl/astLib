#!/usr/bin/env python
#File: cosmologyTest/cosmologyTest.py
#Created: Sat Dec 15 17:23:27 2012
#Last Change: Sat Dec 15 17:25:32 2012
# -*- coding: utf-8 -*-
#
# Tests astCalc routines

from astLib import astCalc
import numpy
import pylab
import IPython
import sys
pylab.matplotlib.interactive(True)

omegaMs = [1.0, 0.2, 0.05]
omegaLs = [0.0, 0.8, 0.0]
styles = ['r-', 'b--', 'g:']
z = numpy.arange(0, 6, 0.1)

# da
pylab.clf()
for m, l, s in zip(omegaMs, omegaLs, styles):
    astCalc.OMEGA_M0 = m
    astCalc.OMEGA_L = l
    label = "$\Omega_{m0}$ = %.2f $\Omega_\Lambda$ = %.2f" % (m, l)
    dH = astCalc.C_LIGHT/astCalc.H0
    plotData = []
    for i in z:
        plotData.append(astCalc.da(i)/dH)
    pylab.plot(z, plotData, s, label=label)
pylab.ylim(0, 0.5)
pylab.xlim(0, 5)
pylab.xlabel("$z$")
pylab.ylabel("$D_A/D_H$")
pylab.legend(loc='lower right')

# dV/dz
pylab.clf()
for m, l, s in zip(omegaMs, omegaLs, styles):
    astCalc.OMEGA_M0 = m
    astCalc.OMEGA_L = l
    label = "$\Omega_{m0}$ = %.2f $\Omega_\Lambda$ = %.2f" % (m, l)
    dH = astCalc.C_LIGHT/astCalc.H0
    plotData = []
    for i in z:
        plotData.append((1.0/dH)**3*astCalc.dVcdz(i))
    pylab.plot(z, plotData, s, label=label)
pylab.ylim(0, 1.1)
pylab.xlim(0, 5)
pylab.xlabel("$z$")
pylab.ylabel("$(1/D_H)^3 dV/dz/d\Omega$")
pylab.legend(loc='upper left')

IPython.embed()
sys.exit()
