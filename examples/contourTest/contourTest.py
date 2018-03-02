#!/usr/bin/env python
#File: contourTest/contourTest.py
#Created: Sat Dec 15 17:04:57 2012
#Last Change: Sat Dec 15 17:05:29 2012
# -*- coding: utf-8 -*-
#
# Test of contour plot

import astropy.io.fits as pyfits
from astLib import *
import pylab
import numpy

TEST_IMAGE = "../../../testingData/testImage2.fits"
TEST_CONTOUR_IMAGE = "../../../testingData/testXRayImage.fits"

img = pyfits.open(TEST_IMAGE)
wcs = astWCS.WCS(TEST_IMAGE)
d = img[0].data

ximg = pyfits.open(TEST_CONTOUR_IMAGE)
xwcs = astWCS.WCS(TEST_CONTOUR_IMAGE)
xd = ximg[0].data

f = astPlots.ImagePlot(d, wcs, axesLabels='decimal')

cLevels = ['log', numpy.median(xd.flatten()), xd.max(), 10]

f.addContourOverlay(xd, xwcs, 'xmm', levels=cLevels, smooth=5.0, color='cyan',
            width=1, highAccuracy=False)
f.draw()
f.save("output_contourTest.png")
