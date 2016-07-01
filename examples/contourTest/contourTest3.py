#!/usr/bin/env python
#File: contourTest/contourTest.py
#Created: Sat Dec 15 17:04:57 2012
#Last Change: Sat Dec 15 17:05:29 2012
# -*- coding: utf-8 -*-
#
# Test of contour plot

import pyfits
from astLib import *
import pylab
import numpy

#TEST_IMAGE = "../../../testingData/testImage2.fits"
TEST_CONTOUR_IMAGE = "../../../../ACTpol/actpol-sourcery/actpol-mfh-sourceryCache/ACT146/ACT-CL_J0000.4+0234.fits"
TEST_IMAGE = TEST_CONTOUR_IMAGE

img = pyfits.open(TEST_IMAGE)
wcs = astWCS.WCS(TEST_IMAGE)
d = img[0].data

ximg = pyfits.open(TEST_CONTOUR_IMAGE)
xwcs = astWCS.WCS(TEST_CONTOUR_IMAGE)
xd = ximg[0].data

pylab.figure(figsize=(40,40))

f = astPlots.ImagePlot(d, wcs, cutLevels = [d.min(), d.max()], axesLabels='decimal', interpolation = "none")

cLevels = numpy.linspace(26.071, 15*26.071, 15)

f.addContourOverlay(xd, xwcs, 'xmm', levels=cLevels, smooth=5.0, color='cyan',
            width=1, highAccuracy=False)
f.draw()
f.save("output_contourTest.png")
