#!/usr/bin/env python
#File: rgbTest/rgbTest.py
#Created: Sat Dec 15 17:29:50 2012
#Last Change: Sat Dec 15 17:35:08 2012
# -*- coding: utf-8 -*-
#
# Example RGB ImagePlot with objects marked and a contour overlay

import numpy
import astropy.io.fits as pyfits
import pylab
from astLib import *

# Load in list of objects to plot - a tab-delimited text file in format id, ra, dec, group
sRAs = []
sDecs = []
sLabels = []
pRAs = []
pDecs = []
pLabels = []
inFile = open("../../../testingData/testData.csv", "r")
lines = inFile.readlines()
for line in lines:
    if line[0] != "#":
        bits = line.split()
        label = bits[0]
        ra = float(bits[1])
        dec = float(bits[2])
        group = int(bits[3])
        if group == 1:
            sLabels.append(label)
            sRAs.append(ra)
            sDecs.append(dec)
        elif group == 2:
            pLabels.append(label)
            pRAs.append(ra)
            pDecs.append(dec)
inFile.close()

# Load the images - these have to be aligned to pixel and same pixel dimensions
rimg = pyfits.open("../../../testingData/testImageR.fits")
gimg = pyfits.open("../../../testingData/testImageG.fits")
bimg = pyfits.open("../../../testingData/testImageB.fits")
wcs = astWCS.WCS("../../../testingData/testImageR.fits")
r = rimg[0].data
g = gimg[0].data
b = bimg[0].data
rCut = [-50, 1000]
gCut = [-50, 1000]
bCut = [-50, 800]

ximg = pyfits.open("../../../testingData/testXRayImage.fits")
xwcs = astWCS.WCS("../../../testingData/testXRayImage.fits")
x = ximg[0].data
cLevels = [2e-05, 2.517e-05, 3.16764e-05, 3.98647e-05, 5.01697e-05,
            6.31385e-05, 7.94597e-05, 0.0001]

# Make the figure - this is a big, high-res version
pylab.figure(figsize=(20,20))
p = astPlots.ImagePlot([r, g, b], wcs, cutLevels = [rCut, gCut, bCut],
    axesLabels="sexagesimal", axesFontSize=26.0, axes = [0.12,0.12,0.8,0.8])

p.addContourOverlay(x, xwcs, 'xmm-contours', levels=cLevels, smooth=5.0)
p.addPlotObjects(sRAs, sDecs, 'spec-members', objLabels=sLabels,
    symbol="circle", color="cyan", width=2.0, size=2.0, objLabelSize=12.0)
p.addPlotObjects(pRAs, pDecs, 'phot-members', objLabels=pLabels, symbol="box",
    color="yellow", width=2.0, size=2.0, objLabelSize=12.0)

# Add a compass because we can
p.addCompass("SW", 10, width=50.0, fontSize=30.0)

p.draw()
p.save("output_rgbTest.png")
