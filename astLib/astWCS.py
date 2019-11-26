"""module for handling World Coordinate Systems (WCS)

(c) 2007-2012 Matt Hilton

(c) 2013-2018 Matt Hilton & Steven Boada

U{http://astlib.sourceforge.net}

This is a higher level interface to some of the routines in PyWCSTools
(distributed with astLib).
PyWCSTools is a simple SWIG wrapping of WCSTools by Jessica Mink
(U{http://tdc-www.harvard.edu/software/wcstools/}). It is intended is to make
this interface complete enough such that direct use of PyWCSTools is
unnecessary.

@var NUMPY_MODE: If True (default), pixel coordinates accepted/returned by
    routines such as L{astWCS.WCS.pix2wcs}, L{astWCS.WCS.wcs2pix} have (0, 0)
    as the origin. Set to False to make these routines accept/return pixel
    coords with (1, 1) as the origin (i.e. to match the FITS convention,
    default behaviour prior to astLib version 0.3.0).
@type NUMPY_MODE: bool

"""

#-----------------------------------------------------------------------------

from astropy.io import fits as pyfits
from PyWCSTools import wcs
import astropy.wcs as apywcs
import numpy
import locale

# if True, -1 from pixel coords to be zero-indexed like numpy. If False, use
# FITS convention.
NUMPY_MODE = True

# Check for the locale bug when decimal separator isn't '.' (atof used in
# libwcs)
lconv = locale.localeconv()
if lconv['decimal_point'] != '.':
    print("WARNING: decimal point separator is not '.' - astWCS coordinate conversions will not work.")
    print("Workaround: after importing any modules that set the locale (e.g. matplotlib) do the following:")
    print("   import locale")
    print("   locale.setlocale(locale.LC_NUMERIC, 'C')")

#-----------------------------------------------------------------------------
class WCS:
    """This class provides methods for accessing information from the World
    Coordinate System (WCS) contained in the header of a FITS image.
    Conversions between pixel and WCS coordinates can also be performed.

    To create a WCS object from a FITS file called "test.fits", simply:

    WCS=astWCS.WCS("test.fits")

    Likewise, to create a WCS object from the pyfits.header of "test.fits":

    img=pyfits.open("test.fits")
    header=img[0].header
    WCS=astWCS.WCS(header, mode = "pyfits")

    """

    def __init__(self, headerSource, extensionName = 0, mode = "image", zapKeywords = [],
                 useAstropyWCS = True):
        """Creates a WCS object using either the information contained in the
        header of the specified .fits image, or from a pyfits.header object.
        Set mode = "pyfits" if the headerSource is a pyfits.header.

        For some images from some archives, particular header keywords such as 
        COMMENT or HISTORY may contain unprintable strings. If you encounter
        this, try setting zapKeywords = ['COMMENT', 'HISTORY'] (for example).
        
        @type headerSource: string or pyfits.header
        @param headerSource: filename of input .fits image, or a pyfits.header
            object
        @type extensionName: int or string
        @param extensionName: name or number of .fits extension in which image
            data is stored
        @type mode: string
        @param mode: set to "image" if headerSource is a .fits file name, or
            set to "astropy" if headerSource is an astropy.io.fits.Header object 
            (setting to this to "pyfits" is equivalent, for backwards 
            compatibility with existing code)
        @type zapKeywords: list
        @param: zapKeywords: keywords to remove from the header before making
            astWCS object.
        @type: useAstropyWCS: bool
        @param: useAstropyWCS: if True, use astropy.wcs to perform WCS 
            coordinate conversions in wcs2pix, pix2wcs (if False, use PyWCSTools)
            
        @note: The meta data provided by headerSource is stored in WCS.header
            as a pyfits.header object.

        """

        self.mode = mode
        self.headerSource = headerSource
        self.extensionName = extensionName

        if self.mode == "image":
            img = pyfits.open(self.headerSource)
            # silentfix below won't deal with unprintable strings
            # so here we optionally remove problematic keywords
            for z in zapKeywords:
                if z in img[self.extensionName].header.keys():
                    for count in range(img[self.extensionName].header.count(z)):
                        img[self.extensionName].header.remove(z)
            img.verify('silentfix') # solves problems with non-standard headers
            self.header = img[self.extensionName].header
            img.close()
        elif self.mode == "astropy" or self.mode == "pyfits":
            for z in zapKeywords:
                if z in self.headerSource.keys():
                    for count in range(self.headerSource.count(z)):
                        self.headerSource.remove(z)
            self.header=headerSource
        
        # If we use astropy.wcs, then we have issues if NAXIS != 2
        if useAstropyWCS == True and self.header['NAXIS'] > 2:
            self.header['NAXIS']=2
        
        # This enables a shim to allow code written for astLib to use astropy.wcs underneath
        self.useAstropyWCS=useAstropyWCS
        if NUMPY_MODE == True:
            self._apywcsOrigin=0
        else:
            self._apywcsOrigin=1

        self.updateFromHeader()


    def copy(self):
        """Copies the WCS object to a new object.

        @rtype: astWCS.WCS object
        @return: WCS object

        """

        # This only sets up a new WCS object, doesn't do a deep copy
        ret = WCS(self.headerSource, self.extensionName, self.mode, 
                  useAstropyWCS = self.useAstropyWCS)

        # This fixes copy bug
        ret.header = self.header.copy()
        ret.updateFromHeader()

        return ret


    def updateFromHeader(self):
        """Updates the WCS object using information from WCS.header. This
        routine should be called whenever changes are made to WCS keywords in
        WCS.header.

        """

        # Updated for pyfits 3.1+
        newHead=pyfits.Header()
        for i in self.header.items():
            if len(str(i[1])) < 70:
                if len(str(i[0])) <= 8:
                    newHead.append((i[0], i[1]))
                else:
                    newHead.append(('HIERARCH '+i[0], i[1]))
        
        # Workaround for ZPN bug when PV2_3 == 0 (as in, e.g., ESO WFI images)
        if "PV2_3" in list(newHead.keys()) and newHead['PV2_3'] == 0 and newHead['CTYPE1'] == 'RA---ZPN':
            newHead["PV2_3"]=1e-15
                
        cardstring = ""
        for card in newHead.cards:
            cardstring = cardstring+str(card)
        
        if self.useAstropyWCS == True:
            self.AWCS = apywcs.WCS(self.header)         # For astropy.wcs shim
        self.WCSStructure = wcs.wcsinit(cardstring)


    def getCentreWCSCoords(self):
        """Returns the RA and dec coordinates (in decimal degrees) at the
        centre of the WCS.

        @rtype: list
        @return: coordinates in decimal degrees in format [RADeg, decDeg]

        """
        full = wcs.wcsfull(self.WCSStructure)
        RADeg = full[0]
        decDeg = full[1]

        return [RADeg, decDeg]


    def getFullSizeSkyDeg(self):
        """Returns the width, height of the image according to the WCS in
        decimal degrees on the sky (i.e., with the projection taken into
        account).

        @rtype: list
        @return: width and height of image in decimal degrees on the sky in
            format [width, height]

        """
        full = wcs.wcsfull(self.WCSStructure)
        width = full[2]
        height = full[3]

        return [width, height]


    def getHalfSizeDeg(self):
        """Returns the half-width, half-height of the image according to the
        WCS in RA and dec degrees.

        @rtype: list
        @return: half-width and half-height of image in R.A., dec. decimal
            degrees in format [half-width, half-height]

        """
        half = wcs.wcssize(self.WCSStructure)
        width = half[2]
        height = half[3]

        return [width, height]


    def getImageMinMaxWCSCoords(self):
        """Returns the minimum, maximum WCS coords defined by the size of the
        parent image (as defined by the NAXIS keywords in the image header).

        @rtype: list
        @return: [minimum R.A., maximum R.A., minimum Dec., maximum Dec.]

        """

        # Get size of parent image this WCS is taken from
        maxX = self.header['NAXIS1']
        maxY = self.header['NAXIS2']
        minX = 1.0
        minY = 1.0

        if NUMPY_MODE == True:
            maxX = maxX-1
            maxY = maxY-1
            minX = minX-1
            minY = minY-1

        bottomLeft = self.pix2wcs(minX, minY)
        topRight = self.pix2wcs(maxX, maxY)

        xCoords = [bottomLeft[0], topRight[0]]
        yCoords = [bottomLeft[1], topRight[1]]
        xCoords.sort()
        yCoords.sort()

        return [xCoords[0], xCoords[1], yCoords[0], yCoords[1]]


    def wcs2pix(self, RADeg, decDeg):
        """Returns the pixel coordinates corresponding to the input WCS
        coordinates (given in decimal degrees). RADeg, decDeg can be single
        floats, or lists or numpy arrays.

        @rtype: list
        @return: pixel coordinates in format [x, y]

        """
        if self.useAstropyWCS == False:
            if type(RADeg) == numpy.ndarray or type(RADeg) == list:
                if type(decDeg) == numpy.ndarray or type(decDeg) == list:
                    pixCoords = []
                    for ra, dec in zip(RADeg, decDeg):
                        pix = wcs.wcs2pix(self.WCSStructure, float(ra), float(dec))
                        # Below handles CEA wraparounds
                        if pix[0] < 1:
                            xTest = ((self.header['CRPIX1'])-(ra-360.0) /
                                self.getXPixelSizeDeg())
                            if xTest >= 1 and xTest < self.header['NAXIS1']:
                                pix[0] = xTest
                        if NUMPY_MODE == True:
                            pix[0] = pix[0]-1
                            pix[1] = pix[1]-1
                        pixCoords.append([pix[0], pix[1]])
            else:
                pixCoords = (wcs.wcs2pix(self.WCSStructure, float(RADeg),
                            float(decDeg)))
                # Below handles CEA wraparounds
                if pixCoords[0] < 1:
                    xTest = ((self.header['CRPIX1'])-(RADeg-360.0) /
                            self.getXPixelSizeDeg())
                    if xTest >= 1 and xTest < self.header['NAXIS1']:
                        pixCoords[0] = xTest
                if NUMPY_MODE == True:
                    pixCoords[0] = pixCoords[0]-1
                    pixCoords[1] = pixCoords[1]-1
                pixCoords = [pixCoords[0], pixCoords[1]]
        
        else:
            # astropy.wcs shim
            pixCoords = self.AWCS.wcs_world2pix(RADeg, decDeg, self._apywcsOrigin)
            pixCoords = numpy.array(pixCoords).transpose().tolist() 

        return pixCoords


    def pix2wcs(self, x, y):
        """Returns the WCS coordinates corresponding to the input pixel
        coordinates.

        @rtype: list
        @return: WCS coordinates in format [RADeg, decDeg]

        """
        if self.useAstropyWCS == False:
            if type(x) == numpy.ndarray or type(x) == list:
                if type(y) == numpy.ndarray or type(y) == list:
                    WCSCoords = []
                    for xc, yc in zip(x, y):
                        if NUMPY_MODE == True:
                            xc += 1
                            yc += 1
                        WCSCoords.append(wcs.pix2wcs(self.WCSStructure, float(xc),
                                    float(yc)))
            else:
                if NUMPY_MODE == True:
                    x += 1
                    y += 1
                WCSCoords = wcs.pix2wcs(self.WCSStructure, float(x), float(y))
        else:
            # astropy.wcs shim
            WCSCoords = self.AWCS.wcs_pix2world(x, y, self._apywcsOrigin)
            WCSCoords = numpy.array(WCSCoords).transpose().tolist() 
            
        return WCSCoords


    def coordsAreInImage(self, RADeg, decDeg):
        """Returns True if the given RA, dec coordinate is within the image
        boundaries.

        @rtype: bool
        @return: True if coordinate within image, False if not.

        """
        
        if self.useAstropyWCS == False:
            pixCoords = wcs.wcs2pix(self.WCSStructure, RADeg, decDeg)
            if pixCoords[0] >= 0 and pixCoords[0] < self.header['NAXIS1'] and \
                pixCoords[1] >= 0 and pixCoords[1] < self.header['NAXIS2']:
                    return True
            else:
                return False
        else:
            # astropy.wcs shim
            import IPython
            print("implement for shim")
            IPython.embed()
            sys.exit()


    def getRotationDeg(self):
        """Returns the rotation angle in degrees around the axis, North through
        East.

        @rtype: float
        @return: rotation angle in degrees

        """
        return self.WCSStructure.rot


    def isFlipped(self):
        """Returns 1 if image is reflected around axis, otherwise returns 0.

        @rtype: int
        @return: 1 if image is flipped, 0 otherwise

        """
        return self.WCSStructure.imflip


    def getPixelSizeDeg(self):
        """Returns the pixel scale of the WCS. This is the average of the x, y
        pixel scales.

        @rtype: float
        @return: pixel size in decimal degrees

        """

        avSize = (abs(self.WCSStructure.xinc)+abs(self.WCSStructure.yinc))/2.0

        return avSize


    def getXPixelSizeDeg(self):
        """Returns the pixel scale along the x-axis of the WCS in degrees.

        @rtype: float
        @return: pixel size in decimal degrees

        """

        avSize = abs(self.WCSStructure.xinc)

        return avSize


    def getYPixelSizeDeg(self):
        """Returns the pixel scale along the y-axis of the WCS in degrees.

        @rtype: float
        @return: pixel size in decimal degrees

        """

        avSize = abs(self.WCSStructure.yinc)

        return avSize


    def getEquinox(self):
        """Returns the equinox of the WCS.

        @rtype: float
        @return: equinox of the WCS

        """
        return self.WCSStructure.equinox


    def getEpoch(self):
        """Returns the epoch of the WCS.

        @rtype: float
        @return: epoch of the WCS

        """
        return self.WCSStructure.epoch


#-----------------------------------------------------------------------------
# Functions for comparing WCS objects
def findWCSOverlap(wcs1, wcs2):
    """Finds the minimum, maximum WCS coords that overlap between wcs1 and
    wcs2. Returns these coordinates, plus the corresponding pixel coordinates
    for each wcs. Useful for clipping overlapping region between two images.

    @rtype: dictionary
    @return: dictionary with keys 'overlapWCS' (min, max RA, dec of overlap
        between wcs1, wcs2) 'wcs1Pix', 'wcs2Pix' (pixel coords in each input
        WCS that correspond to 'overlapWCS' coords)

    """

    mm1 = wcs1.getImageMinMaxWCSCoords()
    mm2 = wcs2.getImageMinMaxWCSCoords()

    overlapWCSCoords = [0.0, 0.0, 0.0, 0.0]

    # Note order swapping below is essential
    # Min RA
    if mm1[0] - mm2[0] <= 0.0:
        overlapWCSCoords[0] = mm2[0]
    else:
        overlapWCSCoords[0] = mm1[0]

    # Max RA
    if mm1[1] - mm2[1] <= 0.0:
        overlapWCSCoords[1] = mm1[1]
    else:
        overlapWCSCoords[1] = mm2[1]

    # Min dec.
    if mm1[2] - mm2[2] <= 0.0:
        overlapWCSCoords[2] = mm2[2]
    else:
        overlapWCSCoords[2] = mm1[2]

    # Max dec.
    if mm1[3] - mm2[3] <= 0.0:
        overlapWCSCoords[3] = mm1[3]
    else:
        overlapWCSCoords[3] = mm2[3]

    # Get corresponding pixel coords
    p1Low = wcs1.wcs2pix(overlapWCSCoords[0], overlapWCSCoords[2])
    p1High = wcs1.wcs2pix(overlapWCSCoords[1], overlapWCSCoords[3])
    p1 = [p1Low[0], p1High[0], p1Low[1], p1High[1]]

    p2Low = wcs2.wcs2pix(overlapWCSCoords[0], overlapWCSCoords[2])
    p2High = wcs2.wcs2pix(overlapWCSCoords[1], overlapWCSCoords[3])
    p2 = [p2Low[0], p2High[0], p2Low[1], p2High[1]]

    return {'overlapWCS': overlapWCSCoords, 'wcs1Pix': p1, 'wcs2Pix': p2}

#-----------------------------------------------------------------------------
