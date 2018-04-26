"""module for coordinate manipulation (conversions, calculations etc.)

(c) 2007-2012 Matt Hilton

(c) 2013-2016 Matt Hilton & Steven Boada

U{http://astlib.sourceforge.net}

"""

import sys
import math
import numpy
from PyWCSTools import wcscon
#import IPython

#-----------------------------------------------------------------------------
def hms2decimal(RAString, delimiter):
    """Converts a delimited string of Hours:Minutes:Seconds format into decimal
    degrees.

    @type RAString: string
    @param RAString: coordinate string in H:M:S format
    @type delimiter: string
    @param delimiter: delimiter character in RAString
    @rtype: float
    @return: coordinate in decimal degrees

    """
    # is it in HH:MM:SS format?
    if delimiter == "":
        RABits = str(RAString).split()
    else:
        RABits = str(RAString).split(delimiter)
    if len(RABits) > 1:
        RAHDecimal = float(RABits[0])
        if len(RABits) > 1:
            RAHDecimal = RAHDecimal+(float(RABits[1])/60.0)
        if len(RABits) > 2:
            RAHDecimal = RAHDecimal+(float(RABits[2])/3600.0)
        RADeg = (RAHDecimal/24.0)*360.0
    else:
        RADeg = float(RAString)

    return RADeg

#-----------------------------------------------------------------------------
def dms2decimal(decString, delimiter):
    """Converts a delimited string of Degrees:Minutes:Seconds format into
    decimal degrees.

    @type decString: string
    @param decString: coordinate string in D:M:S format
    @type delimiter: string
    @param delimiter: delimiter character in decString
    @rtype: float
    @return: coordinate in decimal degrees

    """
    # is it in DD:MM:SS format?
    if delimiter == "":
        decBits = str(decString).split()
    else:
        decBits = str(decString).split(delimiter)
    if len(decBits) > 1:
        decDeg = float(decBits[0])
        if decBits[0].find("-") != -1:
            if len(decBits) > 1:
                decDeg = decDeg-(float(decBits[1])/60.0)
            if len(decBits) > 2:
                decDeg = decDeg-(float(decBits[2])/3600.0)
        else:
            if len(decBits) > 1:
                decDeg = decDeg+(float(decBits[1])/60.0)
            if len(decBits) > 2:
                decDeg = decDeg+(float(decBits[2])/3600.0)
    else:
        decDeg = float(decString)

    return decDeg

#-----------------------------------------------------------------------------
def decimal2hms(RADeg, delimiter):
    """Converts decimal degrees to string in Hours:Minutes:Seconds format with
    user specified delimiter.

    @type RADeg: float
    @param RADeg: coordinate in decimal degrees
    @type delimiter: string
    @param delimiter: delimiter character in returned string
    @rtype: string
    @return: coordinate string in H:M:S format

    """
    hours = (RADeg/360.0)*24
    #if hours < 10 and hours >= 1:
    if 1 <= hours < 10:
        sHours = "0"+str(hours)[0]
    elif hours >= 10:
        sHours = str(hours)[:2]
    elif hours < 1:
        sHours = "00"

    if str(hours).find(".") == -1:
        mins = float(hours)*60.0
    else:
        mins = float(str(hours)[str(hours).index("."):])*60.0
    #if mins<10 and mins>=1:
    if 1 <= mins<10:
        sMins = "0"+str(mins)[:1]
    elif mins >= 10:
        sMins = str(mins)[:2]
    elif mins < 1:
        sMins = "00"

    secs = (hours-(float(sHours)+float(sMins)/60.0))*3600.0
    #if secs < 10 and secs>0.001:
    if 0.001 < secs < 10:
        sSecs = "0"+str(secs)[:str(secs).find(".")+4]
    elif secs < 0.0001:
        sSecs = "00.000"
    else:
        sSecs = str(secs)[:str(secs).find(".")+4]
    if len(sSecs) < 5:
        sSecs = sSecs+"00"	# So all to 3dp

    if float(sSecs) == 60.000:
        sSecs = "00.00"
        sMins = str(int(sMins)+1)
    if int(sMins) == 60:
        sMins = "00"
        sHours = str(int(sHours)+1)

    return sHours+delimiter+sMins+delimiter+sSecs

#------------------------------------------------------------------------------
def decimal2dms(decDeg, delimiter):
    """Converts decimal degrees to string in Degrees:Minutes:Seconds format
    with user specified delimiter.

    @type decDeg: float
    @param decDeg: coordinate in decimal degrees
    @type delimiter: string
    @param delimiter: delimiter character in returned string
    @rtype: string
    @return: coordinate string in D:M:S format

    """
    # Positive
    if decDeg > 0:
        #if decDeg < 10 and decDeg>=1:
        if 1 <= decDeg < 10:
            sDeg = "0"+str(decDeg)[0]
        elif decDeg >= 10:
            sDeg = str(decDeg)[:2]
        elif decDeg < 1:
            sDeg = "00"

        if str(decDeg).find(".") == -1:
            mins = float(decDeg)*60.0
        else:
            mins = float(str(decDeg)[str(decDeg).index("."):])*60
        #if mins<10 and mins>=1:
        if 1 <= mins < 10:
            sMins = "0"+str(mins)[:1]
        elif mins >= 10:
            sMins = str(mins)[:2]
        elif mins < 1:
            sMins = "00"

        secs = (decDeg-(float(sDeg)+float(sMins)/60.0))*3600.0
        #if secs<10 and secs>0:
        if 0 < secs < 10:
            sSecs = "0"+str(secs)[:str(secs).find(".")+3]
        elif secs < 0.001:
            sSecs = "00.00"
        else:
            sSecs = str(secs)[:str(secs).find(".")+3]
        if len(sSecs) < 5:
            sSecs = sSecs+"0"	# So all to 2dp

        if float(sSecs) == 60.00:
            sSecs = "00.00"
            sMins = str(int(sMins)+1)
        if int(sMins) == 60:
            sMins = "00"
            sDeg = str(int(sDeg)+1)

        return "+"+sDeg+delimiter+sMins+delimiter+sSecs

    else:
        #if decDeg>-10 and decDeg<=-1:
        if -10 < decDeg <= -1:
            sDeg = "-0"+str(decDeg)[1]
        elif decDeg <= -10:
            sDeg = str(decDeg)[:3]
        elif decDeg > -1:
            sDeg = "-00"

        if str(decDeg).find(".") == -1:
            mins = float(decDeg)*-60.0
        else:
            mins = float(str(decDeg)[str(decDeg).index("."):])*60
        #if mins<10 and mins>=1:
        if 1 <= mins < 10:
            sMins = "0"+str(mins)[:1]
        elif mins >= 10:
            sMins = str(mins)[:2]
        elif mins < 1:
            sMins = "00"

        secs = (decDeg-(float(sDeg)-float(sMins)/60.0))*3600.0
        #if secs>-10 and secs<0:
        # so don't get minus sign
        if -10 < secs < 0:
            sSecs = "0"+str(secs)[1:str(secs).find(".")+3]
        elif secs > -0.001:
            sSecs = "00.00"
        else:
            sSecs = str(secs)[1:str(secs).find(".")+3]
        if len(sSecs) < 5:
            sSecs = sSecs+"0"	# So all to 2dp

        if float(sSecs) == 60.00:
            sSecs = "00.00"
            sMins = str(int(sMins)+1)
        if int(sMins) == 60:
            sMins = "00"
            sDeg = str(int(sDeg)-1)

        return sDeg+delimiter+sMins+delimiter+sSecs

#-----------------------------------------------------------------------------
def calcAngSepDeg(RADeg1, decDeg1, RADeg2, decDeg2):
    """Calculates the angular separation of two positions on the sky (specified
    in decimal degrees) in decimal degrees. Note that RADeg2, decDeg2 can be numpy
    arrays.

    @type RADeg1: float
    @param RADeg1: R.A. in decimal degrees for position 1
    @type decDeg1: float
    @param decDeg1: dec. in decimal degrees for position 1
    @type RADeg2: float or numpy array
    @param RADeg2: R.A. in decimal degrees for position 2
    @type decDeg2: float or numpy array
    @param decDeg2: dec. in decimal degrees for position 2
    @rtype: float or numpy array, depending upon type of RADeg2, decDeg2
    @return: angular separation in decimal degrees

    """
    
    a=numpy.sin(numpy.radians(decDeg1))*numpy.sin(numpy.radians(decDeg2))+numpy.cos(numpy.radians(decDeg1))*numpy.cos(numpy.radians(decDeg2))*numpy.cos(numpy.radians(RADeg1-RADeg2))
    mask=numpy.greater(a, 1.0)
    if mask.sum() > 0:
        if type(a) == numpy.ndarray:
            a[mask]=1.0
        else:
            a=1.0
    r=numpy.degrees(numpy.arccos(a))
            
    # Above gives nan when RADeg1, decDeg1 == RADeg1, decDeg2
    indexList=numpy.where(numpy.isnan(r) == True)[0]
    tolerance=1e-6
    if len(indexList) > 0:
        for index in indexList:
            if type(r) == numpy.ndarray:
                if type(RADeg2) == numpy.ndarray:
                    if abs(RADeg1 - RADeg2[index]) < tolerance and abs(decDeg1 -decDeg2[index]) < tolerance:
                        r[index]=0.0
                    else:
                        raise Exception("astCoords: calcAngSepDeg - encountered nan not due to equal RADeg, decDeg coords")
                elif type(RADeg1) == numpy.ndarray:
                    if abs(RADeg2 - RADeg1[index]) < tolerance and abs(decDeg2 -decDeg1[index]) < tolerance: 
                        r[index]=0.0
                    else:
                        raise Exception("astCoords: calcAngSepDeg - encountered nan not due to equal RADeg, decDeg coords")
            else:
                r=0.0
        
    return r

#-----------------------------------------------------------------------------
def shiftRADec(ra1, dec1, deltaRA, deltaDec):
    """Computes new right ascension and declination shifted from the original
    by some delta RA and delta DEC. Input position is decimal degrees. Shifts
    (deltaRA, deltaDec) are arcseconds, and output is decimal degrees. Based on
    an IDL routine of the same name.

    @param ra1: float
    @type ra1: R.A. in decimal degrees
    @param dec1: float
    @type dec1: dec. in decimal degrees
    @param deltaRA: float
    @type deltaRA: shift in R.A. in arcseconds
    @param deltaDec: float
    @type deltaDec: shift in dec. in arcseconds
    @rtype: float [newRA, newDec]
    @return: shifted R.A. and dec.

    """

    d2r = math.pi/180.
    as2r = math.pi/648000.

    # Convert everything to radians
    rara1 = ra1*d2r
    dcrad1 = dec1*d2r
    shiftRArad = deltaRA*as2r
    shiftDCrad = deltaDec*as2r

    # Shift!
    deldec2 = 0.0
    sindis = math.sin(shiftRArad / 2.0)
    sindelRA = sindis / math.cos(dcrad1)
    delra = 2.0*math.asin(sindelRA) / d2r

    # Make changes
    ra2 = ra1+delra
    dec2 = dec1 +deltaDec / 3600.0

    return ra2, dec2

#-----------------------------------------------------------------------------
def convertCoords(inputSystem, outputSystem, coordX, coordY, epoch):
    """Converts specified coordinates (given in decimal degrees) between J2000,
    B1950, and Galactic.

    @type inputSystem: string
    @param inputSystem: system of the input coordinates (either "J2000",
        "B1950" or "GALACTIC")
    @type outputSystem: string
    @param outputSystem: system of the returned coordinates (either "J2000",
        "B1950" or "GALACTIC")
    @type coordX: float
    @param coordX: longitude coordinate in decimal degrees, e.g. R. A.
    @type coordY: float
    @param coordY: latitude coordinate in decimal degrees, e.g. dec.
    @type epoch: float
    @param epoch: epoch of the input coordinates
    @rtype: list
    @return: coordinates in decimal degrees in requested output system

    """

    if inputSystem=="J2000" or inputSystem=="B1950" or inputSystem=="GALACTIC":
        if outputSystem=="J2000" or outputSystem=="B1950" or \
            outputSystem=="GALACTIC":

            outCoords=wcscon.wcscon(wcscon.wcscsys(inputSystem),
                wcscon.wcscsys(outputSystem), 0, 0, coordX, coordY, epoch)

            return outCoords

    raise Exception("inputSystem and outputSystem must be 'J2000', 'B1950'"
                    "or 'GALACTIC'")

#-----------------------------------------------------------------------------
def calcRADecSearchBox(RADeg, decDeg, radiusSkyDeg):
    """Calculates minimum and maximum RA, dec coords needed to define a box
    enclosing a circle of radius radiusSkyDeg around the given RADeg, decDeg
    coordinates. Useful for freeform queries of e.g. SDSS, UKIDSS etc.. Uses
    L{calcAngSepDeg}, so has the same limitations.

    @type RADeg: float
    @param RADeg: RA coordinate of centre of search region
    @type decDeg: float
    @param decDeg: dec coordinate of centre of search region
    @type radiusSkyDeg: float
    @param radiusSkyDeg: radius in degrees on the sky used to define search
        region
    @rtype: list
    @return: [RAMin, RAMax, decMin, decMax] - coordinates in decimal degrees
        defining search box

    """

    tolerance = 1e-5  # in degrees on sky
    targetHalfSizeSkyDeg = radiusSkyDeg
    funcCalls = ["calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
               "calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
               "calcAngSepDeg(RADeg, decDeg, RADeg, guess)",
               "calcAngSepDeg(RADeg, decDeg, RADeg, guess)"]
    coords = [RADeg, RADeg, decDeg, decDeg]
    signs = [1.0, -1.0, 1.0, -1.0]
    results = []
    for f, c, sign in zip(funcCalls, coords, signs):
        # Initial guess range
        maxGuess = sign*targetHalfSizeSkyDeg*80.0
        minGuess = sign*targetHalfSizeSkyDeg/80.0
        #guessStep = (maxGuess-minGuess)/10.0
        guesses = numpy.linspace(minGuess+c, maxGuess+c, 1000)
        converged=False
        for i in range(50):
            minSizeDiff = 1e6
            bestGuess = None
            for guess in guesses:
                sizeDiff = abs(eval(f)-targetHalfSizeSkyDeg)
                if sizeDiff < minSizeDiff:
                    minSizeDiff = sizeDiff
                    bestGuess = guess
            if minSizeDiff < tolerance:
                converged=True
                break
            else:
                #print sizeDiff, bestGuess, bestGuess-minGuess, bestGuess-maxGuess
                if bestGuess == None:
                    raise Exception("bestGuess is None")
                guessRange = abs((maxGuess-minGuess))
                maxGuess = bestGuess+guessRange/4.0
                minGuess = bestGuess-guessRange/4.0
                # Stop us from searching the wrong side of the coordinate
                if sign == 1:
                    if minGuess < c:
                        minGuess=c+tolerance
                if sign == -1:
                    if maxGuess > c:
                        maxGuess=c-tolerance
                #guessStep = (maxGuess-minGuess)/20.0
                guesses = numpy.linspace(minGuess, maxGuess, 1000)
        if converged == False:
            raise Exception("calcRADecSearchBox failed to converge")
        results.append(bestGuess)

    RAMax = results[0]
    RAMin = results[1]
    decMax = results[2]
    decMin = results[3]

    # Sanity check
    if (RAMax-RAMin)+(2*tolerance) < 2*targetHalfSizeSkyDeg:
        raise Exception("calcRADecSearchBox failed sanity check")

    return [RAMin, RAMax, decMin, decMax]

