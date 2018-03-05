#!/usr/bin/env python
#File: scaleTest/scaleTest.py
#Created: Sat Dec 15 17:15:17 2012
#Last Change: Sat Dec 15 17:15:56 2012
# -*- coding: utf-8 -*-
#
# Test of image scaling

from astLib import *
import astropy.io.fits as pyfits

TEST_IMAGE = "../../../testingData/testImage1.fits"

img = pyfits.open(TEST_IMAGE)
wcs = astWCS.WCS(TEST_IMAGE)
d = img[0].data

scaled = astImages.scaleImage(d, wcs, 0.55)
astImages.saveFITS("output_scaleTest.fits", scaled['data'], scaled['wcs'])

scaledUp = astImages.scaleImage(d, wcs, 1.43)
astImages.saveFITS("output_scaleTest_up.fits", scaledUp['data'], scaledUp['wcs'])

