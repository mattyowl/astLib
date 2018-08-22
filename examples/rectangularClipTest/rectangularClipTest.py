#!/usr/bin/env python
#File: rectangularClipTest/rectangularClipTest.py
#Created: Sat Dec 15 17:26:52 2012
#Last Change: Sat Dec 15 17:29:07 2012
# -*- coding: utf-8 -*-
#
# Quick script to test out rotation, clipping - rectangles

import astropy.io.fits as pyfits
from astLib import *

# Regular clipping routines test
img = pyfits.open("../../../testingData/testImage3.fits")
wcs = astWCS.WCS("../../../testingData/testImage3.fits")
d = img[0].data

RADeg, decDeg = wcs.getCentreWCSCoords()
RADeg = RADeg+1.2/60.0
decDeg = decDeg-1.0/60.0
widthDeg = 18.0/60.0
heightDeg = 12.3/60.0

clip = astImages.clipRotatedImageSectionWCS(d, wcs, RADeg, decDeg, [widthDeg,
        heightDeg])
clipnr = astImages.clipImageSectionWCS(d, wcs, RADeg, decDeg, [widthDeg,
        heightDeg])
clippix = astImages.clipImageSectionPix(d, 500.0, 500.0, [100.0, 200.0])

astImages.saveFITS("output_rotated.fits", clip['data'], clip['wcs'])
astImages.saveFITS("output_notRotated.fits", clipnr['data'], clipnr['wcs'])
astImages.saveFITS("output_pixelCoords.fits", clippix, None)

# RA, dec coords clipping routine
img = pyfits.open("../../../testingData/testCEAImage.fits")
wcs = astWCS.WCS("../../../testingData/testCEAImage.fits")
d = img[0].data
clip = astImages.clipUsingRADecCoords(img[0].data, wcs, 30.0, 50.0, -55.0,
        -50.0)

astImages.saveFITS("output_RADecClipped.fits", clip['data'], clip['wcs'])
print(clip['wcs'].getImageMinMaxWCSCoords())
