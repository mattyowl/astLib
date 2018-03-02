#!/usr/bin/python
#File: NEDQueryPlot/NEDQueryPlot.py
#Created: Sat Dec 15 17:16:16 2012
#Last Change: Sat Dec 15 17:23:11 2012
# -*- coding: utf-8 -*-
#
# Fetches DSS image at a given position from Skyview, makes a plot listing NED
# matches found in image

#-----------------------------------------------------------------------------
import sys
import os
import math
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve
from astLib import *
import astropy.io.fits as pyfits
import pylab

#-----------------------------------------------------------------------------
def queryNED(RAMin, RAMax, decMin, decMax):
    """Queries NED, returns list of objects. RAMin, RAMax, decMin, decMax, all
    in degrees

    """

    # Just galaxies, groups
    urlretrieve("http://nedwww.ipac.caltech.edu/cgi-bin/nph-allsky?ra_constraint=Between&ra_1=%.6fd&ra_2=%.6fd&dec_constraint=Between&dec_1=%.6fd&dec_2=%.6fd&glon_constraint=Unconstrained&glon_1=&glon_2=&glat_constraint=Unconstrained&glat_1=&glat_2=&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&flux_constraint=Unconstrained&flux_value1=&flux_value2=&flux_unit=Jy&frat_constraint=Unconstrained&ot_include=ANY&in_objtypes1=GGroups&in_objtypes1=Galaxies&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=ascii_tab&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
    % (RAMin, RAMax, decMin, decMax), "result.txt")

    inFile = open("result.txt", "r")
    lines = inFile.readlines()

    dataStarted = False
    labels = []
    names = []
    RAs = []
    decs = []
    sourceTypes = []
    redshifts = []
    for line in lines:
        bits = line.split("\t")
        if bits[0] == "1":
            dataStarted = True
        if dataStarted == True:
            labels.append(bits[0])
            names.append(bits[1])
            RAs.append(float(bits[2]))
            decs.append(float(bits[3]))
            sourceTypes.append(str(bits[4]))
            if bits[6] == '':
                redshifts.append('N/A')
            else:
                redshifts.append(str(bits[6]))

    return ({'labels': labels, 'names': names, 'RAs': RAs, 'decs': decs,
            'sourceTypes': sourceTypes, 'redshifts': redshifts})
    #-------------------------------------------------------------------------
# Main

# Query NED at this position - it's Stephan's quintet
RADeg = 338.98963
decDeg = 33.95991
name = "Stephan's Quintet"

if os.path.exists("fitsImages") == False:
    os.makedirs("fitsImages")
fitsFolder = "fitsImages"

outImageName = name.replace(" ", "_")+".png"
outFITS = name.replace(" ", "_")+".fits"

print("... fetching image from skyview ...")
if os.path.exists(fitsFolder+os.path.sep+outFITS) == True:
    os.remove(fitsFolder+os.path.sep+outFITS)
urlretrieve("http://skyview.gsfc.nasa.gov/cgi-bin/images?Position="+str(RADeg)+","+str(decDeg)+"&Size=0.1&Pixels=1000&Projection=Tan&Grid=J2000&Survey=digitized+sky+survey&Coordinates=J2000&Return=fits",
fitsFolder+os.path.sep+outFITS)

img = pyfits.open(fitsFolder+os.path.sep+outFITS)
wcs = astWCS.WCS(fitsFolder+os.path.sep+outFITS)
imgData = img[0].data
RAMin, RAMax, decMin, decMax=wcs.getImageMinMaxWCSCoords()
objs = queryNED(RAMin, RAMax, decMin, decMax)

fig = pylab.figure(figsize=(10,14))
pylab.figtext(0.5, 0.96, name, fontsize=20, horizontalalignment='center')
pylab.figtext(0.5, 0.94, "R.A, Dec. (J2000) = %s, %s" %
    (astCoords.decimal2hms(RADeg, ":"), astCoords.decimal2dms(decDeg, ":")),
    fontsize=14, horizontalalignment='center')

plot = astPlots.ImagePlot(imgData, wcs, cutLevels=[2000,7000], axes=[0.125,
        0.425, 0.8, 0.5])
plot.addPlotObjects([RADeg], [decDeg], 'actPosition', symbol='cross',
    color='cyan', size=8.0)

pylab.figtext(0.125, 0.35, "NED Matches:", fontweight='bold', fontsize=16)
pylab.figtext(0.125, 0.325, "ID", fontsize=14, horizontalalignment='left',
    verticalalignment='baseline', fontweight='bold')
pylab.figtext(0.16, 0.325, "Name", fontsize=14, horizontalalignment='left',
    verticalalignment='baseline', fontweight='bold')
pylab.figtext(0.42, 0.325, "R.A. (deg.)", fontsize=14,
    horizontalalignment='left', verticalalignment='baseline',
    fontweight='bold')
pylab.figtext(0.57, 0.325, "Dec. (deg.)", fontsize=14,
    horizontalalignment='left', verticalalignment='baseline',
    fontweight='bold')
pylab.figtext(0.72, 0.325, "Type", fontsize=14, horizontalalignment='left',
    verticalalignment='baseline', fontweight='bold')
pylab.figtext(0.8, 0.325, "Redshift", fontsize=14, horizontalalignment='left',
    verticalalignment='baseline', fontweight='bold')
rowSpacing=0.02

if len(objs['RAs']) > 0:
    plot.addPlotObjects(objs['RAs'], objs['decs'], 'nedObjects',
        objLabels=objs['labels'], size=8.0)
    for i in range(len(objs['RAs'])):
        pylab.figtext(0.125, 0.32-(i+1)*rowSpacing, objs['labels'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
        pylab.figtext(0.16, 0.32-(i+1)*rowSpacing, objs['names'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
        pylab.figtext(0.42, 0.32-(i+1)*rowSpacing, "%.6f" % objs['RAs'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
        pylab.figtext(0.57, 0.32-(i+1)*rowSpacing, "%.6f" % objs['decs'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
        pylab.figtext(0.72, 0.32-(i+1)*rowSpacing, objs['sourceTypes'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
        pylab.figtext(0.8, 0.32-(i+1)*rowSpacing, objs['redshifts'][i],
            fontsize=10, horizontalalignment='left',
            verticalalignment='baseline')
else:
    pylab.figtext(0.5, 0.32-rowSpacing, 'No matches found', fontsize=10,
        fontstyle='italic', horizontalalignment='center',
        verticalalignment='baseline')

plot.draw()
plot.save(outImageName)

pylab.close()
img.close()
