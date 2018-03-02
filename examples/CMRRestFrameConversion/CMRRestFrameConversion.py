#!/usr/bin/env python
#File: CMRRestFrameConversion.py
#Created: Sat Dec 15 17:03:04 2012
#Last Change: Sat Dec 15 17:03:50 2012
# -*- coding: utf-8 -*-
#
# Calculates (U-V)z slopes, scatters and intercepts of a given color--magnitude
# relation.
# Follows the procedure in Appendix II of Mei et al. 2009 (ApJ, 690, 42),
# except here we convert
# mags to apparent mags at the distance of the Coma cluster.

from astLib import astSED
from astLib import astStats
from astLib import astCalc
from scipy import stats
from scipy import optimize
from scipy import interpolate
import numpy
import pylab
import os
import sys
import random
import math
import pickle
import string
random.seed()

#-----------------------------------------------------------------------------
# Constants etc.

# number of bootstrap samples, for estimating errors
BOOTSTRAPS = 1000
# number of galaxies, gets * number of models, having same age in each
# metallicity bin
NGAL = 25
FILTER_DIR = "../../../testingData/filters/"

# Map between short filter names on command line and paths, labels etc..
filterMap=[]
filterMap.append({'shortName': 'r625', 'filePath': FILTER_DIR+'F625W_WFC.res',
'plotLabel': 'r625'})
filterMap.append({'shortName': 'i775', 'filePath': FILTER_DIR+'F775W_WFC.res',
'plotLabel': 'i775'})
filterMap.append({'shortName': 'z850', 'filePath': FILTER_DIR+'F850LP_WFC.res',
'plotLabel': 'z850'})
filterMap.append({'shortName': 'U', 'filePath': FILTER_DIR+'U_Johnson.res',
'plotLabel': 'U'})
filterMap.append({'shortName': 'V', 'filePath': FILTER_DIR+'V_Johnson.res',
'plotLabel': 'V'})

# Literature CMR results we want to convert, in their native format
litCMRs = []
litCMRs.append({'name': 'RX J0152.7-1357 (Mei et al. 2009)',
                'redshift': 0.83,
                'colour': 'r625-z850',
                'mag': 'i775',
                'slope': -0.040,
                'slopeErr': 0.017,
                'intercept': 1.93,
                'interceptErr': 0.02,
                'zeroMag': 22.5,
                'scatter': 0.079,
                'scatterErr': 0.008,
                'magType': "AB"})

#-----------------------------------------------------------------------------
def GetCMR(nameFragment, dictList):
    """Finds the CMR dictionary in the litCMRs or results list, by looking for
    nameFragment in name.

    """

    foundCMR = None

    for dict in dictList:
        if nameFragment in dict['name']:
            foundCMR = dict

    return foundCMR

#------------------------------------------------------------------------------
def CalcTransformedCMRWithZM(obsCMR, fitMags, fitCols):
    """Calculates the CMR transformed to the rest frame, using the results of
    the magnitude and colour
    conversion fits. See handwritten notes for the tedious algebra involved.

    """

    # we use the following notation s, zp for slope, zeropoint, append Err for
    # errors
    # CMR for CMR, Mag for fitMags, Col for fitCols
    sCMR = obsCMR['slope']
    zpCMR = obsCMR['intercept']
    zmCMR = obsCMR['zeroMag']

    sMag = fitMags['slope']
    zpMag=fitMags['intercept']

    sCol = fitCols['slope']
    zpCol = fitCols['intercept']

    transformedZeroMag = obsCMR['transformedZeroMag']

    a = sCMR**-1+sMag
    b = sCol/a
    c = zpCMR/sCMR
    d = c-zpMag-zmCMR
    e = d/a
    f = sCol*e
    # this last term is if we want to transform to e.g. match Mei et al.
    g = zpCol+f+(b*transformedZeroMag)

    restCMRSlope = b
    restCMRIntercept = g
    restCMRScatter = obsCMR['scatter']*fitCols['slope']

    return ({'slope': restCMRSlope, 'intercept': restCMRIntercept,
                'scatter': restCMRScatter})

#-----------------------------------------------------------------------------
def BootstrapTransformedCMRErrorsWithZM(obsCMR, fitMags, fitCols):
    """Estimates errors on transformed CMR fit (i.e., into rest frame), by
    assuming the errors on the observed CMR fit, colour transformation fit, and
    mag. transformation fit have Gaussian
    distributions.

    """

    # we use the following notation s, zp for slope, zeropoint, append Err for
    # errors
    # CMR for CMR, Mag for fitMags, Col for fitCols
    sCMR = obsCMR['slope']
    sCMRErr = obsCMR['slopeErr']
    zpCMR = obsCMR['intercept']
    zpCMRErr = obsCMR['interceptErr']

    sMag = fitMags['slope']
    sMagErr = fitMags['slopeError']
    zpMag = fitMags['intercept']
    zpMagErr = fitMags['interceptError']

    sCol = fitCols['slope']
    sColErr = fitCols['slopeError']
    zpCol = fitCols['intercept']
    zpColErr = fitCols['interceptError']

    bsFitResults = []
    for n in range(BOOTSTRAPS):

        bsCMR = {}
        bsMag = {}
        bsCol = {}

        bsCMR['slope'] = random.normalvariate(sCMR, sCMRErr)
        bsCMR['intercept'] = random.normalvariate(zpCMR, zpCMRErr)
        bsCMR['zeroMag'] = obsCMR['zeroMag']
        bsMag['slope'] = random.normalvariate(sMag, sMagErr)
        bsMag['intercept'] = random.normalvariate(zpMag, zpMagErr)
        bsCol['slope'] = random.normalvariate(sCol, sColErr)
        bsCol['intercept'] = random.normalvariate(zpCol, zpColErr)

        bsCMR['scatter'] = random.normalvariate(obsCMR['scatter'],
                        obsCMR['scatterErr'])

        bsCMR['transformedZeroMag'] = obsCMR['transformedZeroMag']

        bsFitResults.append(CalcTransformedCMRWithZM(bsCMR, bsMag, bsCol))

    bsSlopes = []
    bsIntercepts = []
    bsScatters = []
    for bsResult in bsFitResults:
        bsSlopes.append(bsResult['slope'])
        bsIntercepts.append(bsResult['intercept'])
        bsScatters.append(bsResult['scatter'])

    bsSlopes = numpy.array(bsSlopes)
    bsIntercepts = numpy.array(bsIntercepts)
    bsScatters = numpy.array(bsScatters)

    restCMRSlopeErr = numpy.std(bsSlopes)
    restCMRInterceptErr = numpy.std(bsIntercepts)
    restCMRScatterErr = numpy.std(bsScatters)

    return ({'slopeErr': restCMRSlopeErr, 'interceptErr':restCMRInterceptErr,
                'scatterErr': restCMRScatterErr})

#-----------------------------------------------------------------------------
def LoadModels(fileNameList, modelType = "bc03"):
    """Creates a list of stellar population models from the given list of model
    file names.

    """

    models = []
    for f in fileNameList:
        if modelType == "bc03":
            models.append(astSED.BC03Model(f))
        elif modelType == "m05":
            models.append(astSED.M05Model(f))

    return models

#-----------------------------------------------------------------------------
def GetPassbandFileNames(inputColour):
    """Given a mag (e.g. i775) or colour string e.g. r625-z850, lookup the
    appropriate file name(s) in the filterMap, and return the paths in a list.

    """

    bands = inputColour.split("-")
    p = []
    for b in bands:
        p.append(None)
    for row in filterMap:
        for i in range(len(bands)):
            if row['shortName'] == bands[i]:
                p[i]=row['filePath']

    if None in p:
        print("ERROR : couldn't parse colour using filterMap")
        sys.exit()
    else:
        return p

#-----------------------------------------------------------------------------
def LoadPassbands(fileNameList, redshift = None, redshiftPassbands = False):
    """Creates a list of passband objects from the given list of passband file
    names.

    """

    passbands = []
    for f in fileNameList:
        p = astSED.Passband(f)
        if redshiftPassbands == True and redshift != None:
            p.wavelength=p.wavelength*(1.0+redshift)
        passbands.append(p)

    return passbands

#-----------------------------------------------------------------------------
def CalcColourMagTransformation(cmr, restColPassbands, restMagPassband):
    """Calculates the transformation equations needed to convert the given cmr
    into the rest frame at Coma, for the given passbands.

    """

    print((">>> Calculating colour, mag transform for CMR " +
        litCMR['name']+"..."))

    inputColour = cmr['colour']
    inputMag = cmr['mag']
    zCluster = cmr['redshift']

    # Range of formation zs to match Mei et al. 2008
    zfMax = 7.0
    zfMin = 2.0

    observedColPassbandFileNames = GetPassbandFileNames(inputColour)
    observedMagPassbandFileName = GetPassbandFileNames(inputMag)

    observedColLabel = inputColour
    observedMagLabel = inputMag

    # Load stuff
    observedColPassbands = LoadPassbands(observedColPassbandFileNames)
    observedMagPassband = LoadPassbands(observedMagPassbandFileName)[0]

    # Generate galaxy models, we'll hold them all in memory here and use them
    # all in a bit
    print("--> Generating simulated galaxy sample ...")
    restGalaxies = []
    observedGalaxies = []
    for n in range(NGAL):

        print("... n = "+str(n+1)+"/"+str(NGAL)+" ...")

        zfChoice = random.uniform(zfMin, zfMax)
        ageChoice = astCalc.tl(zfChoice)-astCalc.tl(zCluster)

        for i in range(len(models)):

            modelChoice = i
            restGalaxies.append(models[modelChoice].getSED(ageChoice, z=0.02))
            observedGalaxies.append(models[modelChoice].getSED(ageChoice,
                z=zCluster))

    # Fit for colour conversion
    observedColours = []
    restColours = []
    for o, r in zip(observedGalaxies, restGalaxies):
        restColours.append(r.calcColour(restColPassbands[0],
            restColPassbands[1], magType="Vega"))
        observedColours.append(o.calcColour(observedColPassbands[0],
            observedColPassbands[1], magType=cmr['magType']))

    restColours = numpy.array(restColours)
    observedColours = numpy.array(observedColours)

    fitColData = []
    for x, y in zip(observedColours, restColours):
        fitColData.append([x, y])

    fitCols=astStats.OLSFit(fitColData)

    res = restColours-(fitCols['slope']*observedColours+fitCols['intercept'])
    # scatter of residuals, use as fit error
    scatter = astStats.biweightScale(res, 6.0)

    # Fit for mag conversion
    restMinusObservedAppMags = []
    for o, r, obsCol, restCol in zip(observedGalaxies, restGalaxies,
        observedColours, restColours):

        restMinusObservedAppMags.append(r.calcMag(restMagPassband,
            magType="Vega")-o.calcMag(observedMagPassband,
            magType=cmr['magType']))

    restMinusObservedAppMags = numpy.array(restMinusObservedAppMags)

    fitMagData = []
    for x, y in zip(observedColours, restMinusObservedAppMags):
        fitMagData.append([x, y])

    fitMags = astStats.OLSFit(fitMagData)

    trans = {'name': cmr['name'],
            'fitCols': fitCols,
            'fitMags': fitMags,
            'observedColours': observedColours,
            'restColours': restColours,
            'restMinusObservedAppMags': restMinusObservedAppMags}

    return trans

#-----------------------------------------------------------------------------
def GetCSPModel(labelToFind, models, modelLabels):
    """Given a list of models and a matching list of labels, returns the model
    matching the given label.

    """

    foundModel = None
    for m, l in zip(models, modelLabels):
        if l == labelToFind:
            foundModel = m

    return foundModel

#-----------------------------------------------------------------------------
def ApplyCMRTransformation(cmr, trans):
    """Applies the magnitude and colour transformations stored in trans to the
CMR.

    """

    print((">>> Transforming CMR "+cmr['name']+" to Coma rest frame"+
        restColour+" ..."))

    # Check that col, mag transformations and cmrs match up, otherwise
    # something is seriously screwed up
    if cmr['name'] != trans['name']:
        print("ERROR: cmrs and transformation lists not paired!")
        sys.exit()

    fitMags = trans['fitMags']
    fitCols = trans['fitCols']
    transformedCMR = CalcTransformedCMRWithZM(cmr, fitMags, fitCols)
    transformedCMRErrs = BootstrapTransformedCMRErrorsWithZM(cmr, fitMags,
                        fitCols)

    result = {'name': cmr['name'],
            'redshift': cmr['redshift'],
            'colour': restColour,
            'mag': restMag,
            'slope': transformedCMR['slope'],
            'slopeErr': transformedCMRErrs['slopeErr'],
            'intercept': transformedCMR['intercept'],
            'interceptErr': transformedCMRErrs['interceptErr'],
            'zeroMag': cmr['transformedZeroMag'],
            'scatter': transformedCMR['scatter'],
            'scatterErr': transformedCMRErrs['scatterErr']}

    return result

#-----------------------------------------------------------------------------
# Main

# Input parameters
MODEL_TYPE="bc03"

modelFileNames = ["../../../testingData/models/tau0p1Gyr_m42.20",
                "../../../testingData/models/tau0p1Gyr_m52.20",
                "../../../testingData/models/tau0p1Gyr_m62.20",
                "../../../testingData/models/tau0p1Gyr_m72.20"]
models = LoadModels(modelFileNames, modelType = MODEL_TYPE)

# Target colour and mag bands
restColPassbandFileNames = [FILTER_DIR+"U_Johnson.res",
                        FILTER_DIR+"B_Johnson.res"]
restMagPassbandFileName = [FILTER_DIR+"B_Johnson.res"]
restColour = "U-B"    # component of output file name
restMag = "B"
restColLabel = "(U-B)rest"
restMagLabel = "B"
restColPassbands = LoadPassbands(restColPassbandFileNames)
restMagPassband = LoadPassbands(restMagPassbandFileName)[0]

# Evaluate CMR zero point in the rest frame of Coma at intercept of zero
for litCMR in litCMRs:
    litCMR['transformedZeroMag'] = 0.0

# Calculate the colour, mag transformations to take each CMR to the Coma rest frame
transformations = []
for litCMR in litCMRs:
    trans = CalcColourMagTransformation(litCMR, restColPassbands,
            restMagPassband)
    transformations.append(trans)

# Transform the literature CMRs to the rest frame passbands at Coma
results = []
for litCMR, trans in zip(litCMRs, transformations):
    result = ApplyCMRTransformation(litCMR, trans)
    results.append(result)

# Write results to a text file
outFile = open("output_CMRRestConversion.txt", "w")

# Colour, mag transformation
outFile.write("# Color, mag %s rest frame transformation fit coeffs:\n" %
            (restColour))
for t in transformations:
    outFile.write("name = %s\n" % (t['name']))
    outFile.write("# Color transformation:\n")
    outFile.write("slope = %.5f\n" % (t['fitCols']['slope']))
    outFile.write("slopeError = %.5f\n" % (t['fitCols']['slopeError']))
    outFile.write("intercept = %.5f\n" % (t['fitCols']['intercept']))
    outFile.write("interceptError = %.5f\n" % (t['fitCols']['interceptError']))
    outFile.write("# Mag transformation:\n")
    outFile.write("slope = %.5f\n" % (t['fitMags']['slope']))
    outFile.write("slopeError = %.5f\n" % (t['fitMags']['slopeError']))
    outFile.write("intercept = %.5f\n" % (t['fitMags']['intercept']))
    outFile.write("interceptError = %.5f\n" % (t['fitMags']['interceptError']))

# Transformed CMR
outFile.write("# Transformed CMR:\n")
keyOrder = ["name", "redshift", "colour", "mag", "zeroMag", "slope",
            "slopeErr", "intercept", "interceptErr", "scatter", "scatterErr"]
for r in results:
    for k in keyOrder:
        for key in list(r.keys()):
            if str(key) == k:
                if type(r[key]) == str:
                    outFile.write("%s = %s\n" % (key, r[key]))
                else:
                    outFile.write("%s = %.3f\n" % (key, float(r[key])))
outFile.close()
