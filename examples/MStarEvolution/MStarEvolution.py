#!/usr/bin/env python
#File: MStarEvolution/MStarEvolution.py
#Created: Sat Dec 15 17:06:36 2012
#Last Change: Sat Dec 15 17:08:15 2012
# -*- coding: utf-8 -*-
#
# Calculates evolution of the characteristic magnitude M* in the galaxy
# luminosity function, in the K band, normalised using De Propris et al. 1999
# (not fitted), and overplots data

from astLib import astSED
import numpy
import pylab

bc03 = astSED.BC03Model("../../../testingData/models/tau0p1Gyr_m62.1")
K = astSED.Passband("../../../testingData/filters/K_2MASS.res")

DP99Mags = [14.84, 15.16, 15.64, 15.74, 16.50, 16.38, 16.85, 17.57, 17.51,
            18.05]
DP99zs = [0.15,  0.20,  0.25,  0.32,  0.40,  0.46,  0.54,  0.61,  0.79,  0.9]
pylab.plot(DP99zs, DP99Mags, 'ro', label='DP1999')

magTrack2 = bc03.getMagEvolution(K, 15.74, 0.32, 2.0, zStepSize=0.1,
            onePlusZSteps=True)
magTrack5 = bc03.getMagEvolution(K, 15.74, 0.32, 5.0, zStepSize=0.1,
            onePlusZSteps=True)

pylab.plot(magTrack2['z'], magTrack2['mag'], 'b--', label='zf=2.0')
pylab.plot(magTrack5['z'], magTrack5['mag'], 'r--', label='zf=5.0')

pylab.xlim(0, 2)
pylab.ylim(13, 22)
pylab.xlabel('redshift')
pylab.ylabel('Ks (Vega)')
pylab.legend()
pylab.savefig("magEvo.png")

