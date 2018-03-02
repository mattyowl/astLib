#!/usr/bin/env python
#File: compassTest/compassTest.py
#Created: Sat Dec 15 17:29:28 2012
#Last Change: Sat Dec 15 17:29:42 2012
# -*- coding: utf-8 -*-
#
# Test of contour plot

import astropy.io.fits as pyfits
from astLib import *
import pylab

TEST_IMAGE = "../../../testingData/testImage1.fits"

img = pyfits.open(TEST_IMAGE)
wcs = astWCS.WCS(TEST_IMAGE)
d = img[0].data

f = astPlots.ImagePlot(d, wcs, axes = [0, 0, 1, 1])
f.addCompass("SE", 60.0, color = 'red')

f.draw()
f.save("output_compassTest.png")
