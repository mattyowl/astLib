#!/usr/bin/env python
#File: simpleRGB/simpleRGB.py
#Created: Sat Dec 15 17:25:44 2012
#Last Change: Sat Dec 15 17:26:38 2012
# -*- coding: utf-8 -*-
#
# Simple RGB ImagePlot example

import numpy
import astropy.io.fits as pyfits
import pylab
from astLib import *

# Load the images - these have to be aligned to pixel and same pixel
# dimensions
r = pyfits.getdata("../../../testingData/stephanDSS2IR.fits")
g = pyfits.getdata("../../../testingData/stephanDSS2Red.fits")
b = pyfits.getdata("../../../testingData/stephanDSS2Blue.fits")
wcs = astWCS.WCS("../../../testingData/stephanDSS2Blue.fits")
rCut = [r.min(), r.max()]
gCut = [g.min(), g.max()]
bCut = [b.min(), b.max()]

# Make the figure
p = astPlots.ImagePlot([r, g, b], wcs, axes = [0.12, 0.12, 0.8, 0.8], cutLevels = [rCut, gCut, bCut],
    title="Stephan's Quintet")
p.save("output_simpleRGB.png")
