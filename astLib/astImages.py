"""module for simple .fits image tasks (rotation, clipping out sections, making .pngs etc.)

(c) 2007-2018 Matt Hilton 

U{http://astlib.sourceforge.net}

Some routines in this module will fail if, e.g., asked to clip a section from a .fits image at a
position not found within the image (as determined using the WCS). Where this occurs, the function
will return None. An error message will be printed to the console when this happens if
astImages.REPORT_ERRORS=True (the default). Testing if an astImages function returns None can be
used to handle errors in scripts. 

"""

REPORT_ERRORS=True

import os
import sys
import math
from astLib import astWCS
from astropy.io import fits as pyfits    
try:
    from scipy import ndimage
    from scipy import interpolate
except ImportError:
    print("WARNING: astImages: failed to import scipy.ndimage - some functions will not work.")
import numpy as np
try:
    import matplotlib
    from matplotlib import pylab
    matplotlib.interactive(False)
except ImportError:
    print("WARNING: astImages: failed to import matplotlib - some functions will not work.")

#---------------------------------------------------------------------------------------------------
def clipImageSectionWCS(imageData, imageWCS, RADeg, decDeg, clipSizeDeg, returnWCS = True):
    """Clips a square or rectangular section from an image array at the given celestial coordinates. 
    An updated WCS for the clipped section is optionally returned, as well as the x, y pixel 
    coordinates in the original image corresponding to the clipped section.
    
    Note that the clip size is specified in degrees on the sky. For projections that have varying
    real pixel scale across the map (e.g. CEA), use L{clipUsingRADecCoords} instead.

    Similarly, this routine will not work for a WCS that has polynomial distortion coefficients 
    in the header (e.g., CTYPE1 = 'RA---TAN-SIP' etc.) - again L{clipUsingRADecCoords} can be used
    in such cases.
    
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS
    @param imageWCS: astWCS.WCS object
    @type RADeg: float
    @param RADeg: coordinate in decimal degrees
    @type decDeg: float
    @param decDeg: coordinate in decimal degrees
    @type clipSizeDeg: float or list in format [widthDeg, heightDeg]
    @param clipSizeDeg: if float, size of square clipped section in decimal degrees; if list,
    size of clipped section in degrees in x, y axes of image respectively
    @type returnWCS: bool
    @param returnWCS: if True, return an updated WCS for the clipped section
    @rtype: dictionary
    @return: clipped image section (np array), updated astWCS WCS object for
    clipped image section, and coordinates of clipped section in imageData in format 
    {'data', 'wcs', 'clippedSection'}.
        
    """	
    
    imHeight=imageData.shape[0]
    imWidth=imageData.shape[1]
    xImScale=imageWCS.getXPixelSizeDeg()
    yImScale=imageWCS.getYPixelSizeDeg()
    
    if type(clipSizeDeg) == float:
        xHalfClipSizeDeg=clipSizeDeg/2.0
        yHalfClipSizeDeg=xHalfClipSizeDeg
    elif type(clipSizeDeg) == list or type(clipSizeDeg) == tuple:
        xHalfClipSizeDeg=clipSizeDeg[0]/2.0
        yHalfClipSizeDeg=clipSizeDeg[1]/2.0
    else:
        raise Exception("did not understand clipSizeDeg: should be float, or [widthDeg, heightDeg]")
    
    xHalfSizePix=xHalfClipSizeDeg/xImScale
    yHalfSizePix=yHalfClipSizeDeg/yImScale    
    
    cPixCoords=imageWCS.wcs2pix(RADeg, decDeg)
    
    cTopLeft=[cPixCoords[0]+xHalfSizePix, cPixCoords[1]+yHalfSizePix]
    cBottomRight=[cPixCoords[0]-xHalfSizePix, cPixCoords[1]-yHalfSizePix]
        
    X=[int(round(cTopLeft[0])),int(round(cBottomRight[0]))]
    Y=[int(round(cTopLeft[1])),int(round(cBottomRight[1]))]
    
    X.sort()
    Y.sort()
    
    if X[0] < 0:
        X[0]=0
    if X[1] > imWidth:
        X[1]=imWidth
    if Y[0] < 0:
        Y[0]=0
    if Y[1] > imHeight:
        Y[1]=imHeight
    
    clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]

    # Update WCS
    if returnWCS == True:
        try:
            oldCRPIX1=imageWCS.header['CRPIX1']
            oldCRPIX2=imageWCS.header['CRPIX2']
            clippedWCS=imageWCS.copy()
            clippedWCS.header['NAXIS1']=clippedData.shape[1]
            clippedWCS.header['NAXIS2']=clippedData.shape[0]
            clippedWCS.header['CRPIX1']=oldCRPIX1-X[0]
            clippedWCS.header['CRPIX2']=oldCRPIX2-Y[0]
            clippedWCS.updateFromHeader()
            
        except KeyError:
            
            if REPORT_ERRORS == True:
                
                print("WARNING: astImages.clipImageSectionWCS() : no CRPIX1, CRPIX2 keywords found - not updating clipped image WCS.")
                
                clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]
                clippedWCS=imageWCS.copy()
    else:
        clippedWCS=None
    
    return {'data': clippedData, 'wcs': clippedWCS, 'clippedSection': [X[0], X[1], Y[0], Y[1]]}
    
#---------------------------------------------------------------------------------------------------
def clipImageSectionPix(imageData, XCoord, YCoord, clipSizePix):
    """Clips a square or rectangular section from an image array at the given pixel coordinates.
    
    @type imageData: np array
    @param imageData: image data array
    @type XCoord: float
    @param XCoord: coordinate in pixels
    @type YCoord: float
    @param YCoord: coordinate in pixels
    @type clipSizePix: float or list in format [widthPix, heightPix]
    @param clipSizePix: if float, size of square clipped section in pixels; if list,
    size of clipped section in pixels in x, y axes of output image respectively
    @rtype: np array
    @return: clipped image section
    
    """		
    
    imHeight=imageData.shape[0]
    imWidth=imageData.shape[1]
    
    if type(clipSizePix) == float or type(clipSizePix) == int:
        xHalfClipSizePix=int(round(clipSizePix/2.0))
        yHalfClipSizePix=xHalfClipSizePix
    elif type(clipSizePix) == list or type(clipSizePix) == tuple:
        xHalfClipSizePix=int(round(clipSizePix[0]/2.0))
        yHalfClipSizePix=int(round(clipSizePix[1]/2.0))
    else:
        raise Exception("did not understand clipSizePix: should be float, or [widthPix, heightPix]")
       
    cTopLeft=[XCoord+xHalfClipSizePix, YCoord+yHalfClipSizePix]
    cBottomRight=[XCoord-xHalfClipSizePix, YCoord-yHalfClipSizePix]
    
    X=[int(round(cTopLeft[0])),int(round(cBottomRight[0]))]
    Y=[int(round(cTopLeft[1])),int(round(cBottomRight[1]))]
    
    X.sort()
    Y.sort()
    
    if X[0] < 0:
        X[0]=0
    if X[1] > imWidth:
        X[1]=imWidth
    if Y[0] < 0:
        Y[0]=0
    if Y[1] > imHeight:
        Y[1]=imHeight		
        
    return imageData[Y[0]:Y[1],X[0]:X[1]]
    
#---------------------------------------------------------------------------------------------------
def clipRotatedImageSectionWCS(imageData, imageWCS, RADeg, decDeg, clipSizeDeg, returnWCS = True):
    """Clips a square or rectangular section from an image array at the given celestial coordinates. 
    The resulting clip is rotated and/or flipped such that North is at the top, and East appears at
    the left. An updated WCS for the clipped section is also returned. Note that the alignment
    of the rotated WCS is currently not perfect - however, it is probably good enough in most
    cases for use with L{ImagePlot} for plotting purposes.
    
    Note that the clip size is specified in degrees on the sky. For projections that have varying
    real pixel scale across the map (e.g. CEA), use L{clipUsingRADecCoords} instead.
    
    Similarly, this routine will not work for a WCS that has polynomial distortion coefficients 
    in the header (e.g., CTYPE1 = 'RA---TAN-SIP' etc.) - again L{clipUsingRADecCoords} can be used
    in such cases.
    
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS
    @param imageWCS: astWCS.WCS object
    @type RADeg: float
    @param RADeg: coordinate in decimal degrees
    @type decDeg: float
    @param decDeg: coordinate in decimal degrees
    @type clipSizeDeg: float
    @param clipSizeDeg: if float, size of square clipped section in decimal degrees; if list,
    size of clipped section in degrees in RA, dec. axes of output rotated image respectively
    @type returnWCS: bool
    @param returnWCS: if True, return an updated WCS for the clipped section
    @rtype: dictionary
    @return: clipped image section (np array), updated astWCS WCS object for
    clipped image section, in format {'data', 'wcs'}.
    
    @note: Returns 'None' if the requested position is not found within the image. If the image
    WCS does not have keywords of the form CD1_1 etc., the output WCS will not be rotated.
    
    """
        
    halfImageSize=imageWCS.getHalfSizeDeg()
    imageCentre=imageWCS.getCentreWCSCoords()
    imScale=imageWCS.getPixelSizeDeg()

    if type(clipSizeDeg) == float:
        xHalfClipSizeDeg=clipSizeDeg/2.0
        yHalfClipSizeDeg=xHalfClipSizeDeg
    elif type(clipSizeDeg) == list or type(clipSizeDeg) == tuple:
        xHalfClipSizeDeg=clipSizeDeg[0]/2.0
        yHalfClipSizeDeg=clipSizeDeg[1]/2.0
    else:
        raise Exception("did not understand clipSizeDeg: should be float, or [widthDeg, heightDeg]")
    
    diagonalHalfSizeDeg=math.sqrt((xHalfClipSizeDeg*xHalfClipSizeDeg) \
        +(yHalfClipSizeDeg*yHalfClipSizeDeg))
    
    diagonalHalfSizePix=diagonalHalfSizeDeg/imScale
        
    if RADeg>imageCentre[0]-halfImageSize[0] and RADeg<imageCentre[0]+halfImageSize[0] \
        and decDeg>imageCentre[1]-halfImageSize[1] and decDeg<imageCentre[1]+halfImageSize[1]:
        
        imageDiagonalClip=clipImageSectionWCS(imageData, imageWCS, RADeg,
                        decDeg, diagonalHalfSizeDeg*2.0)
        diagonalClip=imageDiagonalClip['data']
        diagonalWCS=imageDiagonalClip['wcs']
        
        rotDeg=diagonalWCS.getRotationDeg()
        imageRotated=ndimage.rotate(diagonalClip, rotDeg)
        if diagonalWCS.isFlipped() == 1:
            imageRotated=pylab.fliplr(imageRotated)
        
        # Handle WCS rotation
        rotatedWCS=diagonalWCS.copy()
        rotRadians=math.radians(rotDeg)

        if returnWCS == True:
            try:
                
                CD11=rotatedWCS.header['CD1_1']
                CD21=rotatedWCS.header['CD2_1']
                CD12=rotatedWCS.header['CD1_2']
                CD22=rotatedWCS.header['CD2_2']
                if rotatedWCS.isFlipped() == 1:
                    CD11=CD11*-1
                    CD12=CD12*-1
                CDMatrix=np.array([[CD11, CD12], [CD21, CD22]], dtype=np.float64)

                rotRadians=rotRadians
                rot11=math.cos(rotRadians)
                rot12=math.sin(rotRadians)
                rot21=-math.sin(rotRadians)
                rot22=math.cos(rotRadians)
                rotMatrix=np.array([[rot11, rot12], [rot21, rot22]], dtype=np.float64)
                newCDMatrix=np.dot(rotMatrix, CDMatrix)

                P1=diagonalWCS.header['CRPIX1']
                P2=diagonalWCS.header['CRPIX2']
                V1=diagonalWCS.header['CRVAL1']
                V2=diagonalWCS.header['CRVAL2']
                
                PMatrix=np.zeros((2,), dtype = np.float64)
                PMatrix[0]=P1
                PMatrix[1]=P2
                
                # BELOW IS HOW TO WORK OUT THE NEW REF PIXEL
                CMatrix=np.array([imageRotated.shape[1]/2.0, imageRotated.shape[0]/2.0])
                centreCoords=diagonalWCS.getCentreWCSCoords()
                alphaRad=math.radians(centreCoords[0])
                deltaRad=math.radians(centreCoords[1])
                thetaRad=math.asin(math.sin(deltaRad)*math.sin(math.radians(V2)) + \
                                math.cos(deltaRad)*math.cos(math.radians(V2))*math.cos(alphaRad-math.radians(V1)))
                phiRad=math.atan2(-math.cos(deltaRad)*math.sin(alphaRad-math.radians(V1)), \
                                math.sin(deltaRad)*math.cos(math.radians(V2)) - \
                                math.cos(deltaRad)*math.sin(math.radians(V2))*math.cos(alphaRad-math.radians(V1))) + \
                                math.pi
                RTheta=(180.0/math.pi)*(1.0/math.tan(thetaRad))
                
                xy=np.zeros((2,), dtype=np.float64)
                xy[0]=RTheta*math.sin(phiRad)
                xy[1]=-RTheta*math.cos(phiRad)
                newPMatrix=CMatrix - np.dot(np.linalg.inv(newCDMatrix), xy)
                
                # But there's a small offset to CRPIX due to the rotatedImage being rounded to an integer
                # number of pixels (not sure this helps much)
                #d=np.dot(rotMatrix, [diagonalClip.shape[1], diagonalClip.shape[0]])
                #offset=abs(d)-np.array(imageRotated.shape)
                
                rotatedWCS.header['NAXIS1']=imageRotated.shape[1]
                rotatedWCS.header['NAXIS2']=imageRotated.shape[0]
                rotatedWCS.header['CRPIX1']=newPMatrix[0]
                rotatedWCS.header['CRPIX2']=newPMatrix[1]
                rotatedWCS.header['CRVAL1']=V1
                rotatedWCS.header['CRVAL2']=V2
                rotatedWCS.header['CD1_1']=newCDMatrix[0][0]
                rotatedWCS.header['CD2_1']=newCDMatrix[1][0]
                rotatedWCS.header['CD1_2']=newCDMatrix[0][1]
                rotatedWCS.header['CD2_2']=newCDMatrix[1][1]
                rotatedWCS.updateFromHeader()
                                
            except KeyError:
                
                if REPORT_ERRORS == True:
                    print("WARNING: astImages.clipRotatedImageSectionWCS() : no CDi_j keywords found - not rotating WCS.")
                    
                imageRotated=diagonalClip
                rotatedWCS=diagonalWCS
            
        imageRotatedClip=clipImageSectionWCS(imageRotated, rotatedWCS, RADeg, decDeg, clipSizeDeg)
        
        if returnWCS == True:
            return {'data': imageRotatedClip['data'], 'wcs': imageRotatedClip['wcs']}
        else:
            return {'data': imageRotatedClip['data'], 'wcs': None}
        
    else:
        
        if REPORT_ERRORS==True:
            print("""ERROR: astImages.clipRotatedImageSectionWCS() : 
            RADeg, decDeg are not within imageData.""")
        
        return None

#---------------------------------------------------------------------------------------------------
def clipUsingRADecCoords(imageData, imageWCS, RAMin, RAMax, decMin, decMax, returnWCS = True):
    """Clips a section from an image array at the pixel coordinates corresponding to the given
    celestial coordinates.
    
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS
    @param imageWCS: astWCS.WCS object
    @type RAMin: float
    @param RAMin: minimum RA coordinate in decimal degrees
    @type RAMax: float
    @param RAMax: maximum RA coordinate in decimal degrees
    @type decMin: float
    @param decMin: minimum dec coordinate in decimal degrees
    @type decMax: float
    @param decMax: maximum dec coordinate in decimal degrees
    @type returnWCS: bool
    @param returnWCS: if True, return an updated WCS for the clipped section
    @rtype: dictionary
    @return: clipped image section (np array), updated astWCS WCS object for
    clipped image section, and corresponding pixel coordinates in imageData in format 
    {'data', 'wcs', 'clippedSection'}.
    
    @note: Returns 'None' if the requested position is not found within the image.
    
    """
    
    imHeight=imageData.shape[0]
    imWidth=imageData.shape[1]
    
    # Fixed for TPV headers
    xMin, yMin=imageWCS.wcs2pix(RAMax, decMin)
    xMax, yMax=imageWCS.wcs2pix(RAMin, decMax)
    X=[xMin, xMax]
    X.sort()
    Y=[yMin, yMax]
    Y.sort()
    X=[int(np.floor(X[0])), int(np.ceil(X[1]))]
    Y=[int(np.floor(Y[0])), int(np.ceil(Y[1]))]
    
    if X[0] < 0:
        X[0]=0
    if X[1] > imWidth:
        X[1]=imWidth
    if Y[0] < 0:
        Y[0]=0
    if Y[1] > imHeight:
        Y[1]=imHeight   
    
    clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]

    # Update WCS
    if returnWCS == True:
        try:
            oldCRPIX1=imageWCS.header['CRPIX1']
            oldCRPIX2=imageWCS.header['CRPIX2']
            clippedWCS=imageWCS.copy()
            clippedWCS.header['NAXIS1']=clippedData.shape[1]
            clippedWCS.header['NAXIS2']=clippedData.shape[0]
            clippedWCS.header['CRPIX1']=oldCRPIX1-X[0]
            clippedWCS.header['CRPIX2']=oldCRPIX2-Y[0]
            clippedWCS.updateFromHeader()
            
        except KeyError:
            
            if REPORT_ERRORS == True:
                
                print("WARNING: astImages.clipUsingRADecCoords() : no CRPIX1, CRPIX2 keywords found - not updating clipped image WCS.")
                
                clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]
                clippedWCS=imageWCS.copy()
    else:
        clippedWCS=None
    
    return {'data': clippedData, 'wcs': clippedWCS, 'clippedSection': [X[0], X[1], Y[0], Y[1]]}
    
#---------------------------------------------------------------------------------------------------
def scaleImage(imageData, imageWCS, scaleFactor):
    """Scales image array and WCS by the given scale factor.
    
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS
    @param imageWCS: astWCS.WCS object
    @type scaleFactor: float or list or tuple
    @param scaleFactor: factor to resize image by - if tuple or list, in format 
        [x scale factor, y scale factor]
    @rtype: dictionary
    @return: image data (np array), updated astWCS WCS object for image, in format {'data', 'wcs'}.
    
    """

    if type(scaleFactor) == int or type(scaleFactor) == float:
        scaleFactor=[float(scaleFactor), float(scaleFactor)]    
    scaledData=ndimage.zoom(imageData, scaleFactor)
    
    # Changed below because ndimage.zoom now uses round instead of int (since scipy 0.13.0)
    # NOTE: np axes order flips order compared to scaleFactor
    trueScaleFactor=np.array(scaledData.shape, dtype = float) / np.array(imageData.shape, dtype = float)
    offset=0.
    
    # Rescale WCS
    try:
        oldCRPIX1=imageWCS.header['CRPIX1']
        oldCRPIX2=imageWCS.header['CRPIX2']
        CD11=imageWCS.header['CD1_1']
        CD21=imageWCS.header['CD2_1']
        CD12=imageWCS.header['CD1_2']
        CD22=imageWCS.header['CD2_2'] 
    except KeyError:
        # Try the older FITS header format
        try:
            oldCRPIX1=imageWCS.header['CRPIX1']
            oldCRPIX2=imageWCS.header['CRPIX2']
            CD11=imageWCS.header['CDELT1']
            CD21=0
            CD12=0
            CD22=imageWCS.header['CDELT2']
        except KeyError:
            if REPORT_ERRORS == True:
                print("WARNING: astImages.rescaleImage() : no CDij or CDELT keywords found - not updating WCS.")
            scaledWCS=imageWCS.copy()
            return {'data': scaledData, 'wcs': scaledWCS}

    CDMatrix=np.array([[CD11, CD12], [CD21, CD22]], dtype=np.float64)
    scaleFactorMatrix=np.array([[1.0/trueScaleFactor[1], 0], [0, 1.0/trueScaleFactor[0]]])
    scaleFactorMatrix=np.array([[1.0/trueScaleFactor[1], 0], [0, 1.0/trueScaleFactor[0]]])
    scaledCDMatrix=np.dot(scaleFactorMatrix, CDMatrix)

    scaledWCS=imageWCS.copy()
    scaledWCS.header['NAXIS1']=scaledData.shape[1]
    scaledWCS.header['NAXIS2']=scaledData.shape[0]
    scaledWCS.header['CRPIX1']=oldCRPIX1*trueScaleFactor[1]
    scaledWCS.header['CRPIX2']=oldCRPIX2*trueScaleFactor[0]
    scaledWCS.header['CD1_1']=scaledCDMatrix[0][0]
    scaledWCS.header['CD2_1']=scaledCDMatrix[1][0]
    scaledWCS.header['CD1_2']=scaledCDMatrix[0][1]
    scaledWCS.header['CD2_2']=scaledCDMatrix[1][1]
    scaledWCS.updateFromHeader()
    
    return {'data': scaledData, 'wcs': scaledWCS}
    
#---------------------------------------------------------------------------------------------------
def intensityCutImage(imageData, cutLevels):
    """Creates a matplotlib.pylab plot of an image array with the specified cuts in intensity
    applied. This routine is used by L{saveBitmap} and L{saveContourOverlayBitmap}, which both
    produce output as .png, .jpg, etc. images.
    
    @type imageData: np array
    @param imageData: image data array
    @type cutLevels: list
    @param cutLevels: sets the image scaling - available options:
        - pixel values: cutLevels=[low value, high value].
        - histogram equalisation: cutLevels=["histEq", number of bins ( e.g. 1024)]
        - relative: cutLevels=["relative", cut per cent level (e.g. 99.5)]
        - smart: cutLevels=["smart", cut per cent level (e.g. 99.5)]
    ["smart", 99.5] seems to provide good scaling over a range of different images.
    @rtype: dictionary
    @return: image section (np.array), matplotlib image normalisation (matplotlib.colors.Normalize), in the format {'image', 'norm'}.
    
    @note: If cutLevels[0] == "histEq", then only {'image'} is returned.
    
    """
    
    oImWidth=imageData.shape[1]
    oImHeight=imageData.shape[0]
                    
    # Optional histogram equalisation
    if cutLevels[0]=="histEq":
        
        imageData=histEq(imageData, cutLevels[1])
        anorm=pylab.Normalize(imageData.min(), imageData.max())
        
    elif cutLevels[0]=="relative":
        
        # this turns image data into 1D array then sorts
        sorted=np.sort(np.ravel(imageData))	
        maxValue=sorted.max()
        minValue=sorted.min()
        
        # want to discard the top and bottom specified
        topCutIndex=len(sorted-1) \
            -int(math.floor(float((100.0-cutLevels[1])/100.0)*len(sorted-1)))
        bottomCutIndex=int(math.ceil(float((100.0-cutLevels[1])/100.0)*len(sorted-1)))
        topCut=sorted[topCutIndex]
        bottomCut=sorted[bottomCutIndex]
        anorm=pylab.Normalize(bottomCut, topCut)
        
    elif cutLevels[0]=="smart":
        
        # this turns image data into 1Darray then sorts
        sorted=np.sort(np.ravel(imageData))	
        maxValue=sorted.max()
        minValue=sorted.min()
        numBins=10000 		# 0.01 per cent accuracy
        binWidth=(maxValue-minValue)/float(numBins)
        histogram=ndimage.histogram(sorted, minValue, maxValue, numBins)
        
        # Find the bin with the most pixels in it, set that as our minimum
        # Then search through the bins until we get to a bin with more/or the same number of
        # pixels in it than the previous one.
        # We take that to be the maximum.
        # This means that we avoid the traps of big, bright, saturated stars that cause
        # problems for relative scaling
        backgroundValue=histogram.max()
        foundBackgroundBin=False
        foundTopBin=False
        lastBin=-10000					
        for i in range(len(histogram)):
            
            if histogram[i]>=lastBin and foundBackgroundBin==True:
                
                # Added a fudge here to stop us picking for top bin a bin within 
                # 10 percent of the background pixel value
                if (minValue+(binWidth*i))>bottomBinValue*1.1:
                    topBinValue=minValue+(binWidth*i)
                    foundTopBin=True
                    break
            
            if histogram[i]==backgroundValue and foundBackgroundBin==False:
                bottomBinValue=minValue+(binWidth*i)
                foundBackgroundBin=True

            lastBin=histogram[i]
        
        if foundTopBin==False:
            topBinValue=maxValue
         
        #Now we apply relative scaling to this
        smartClipped=np.clip(sorted, bottomBinValue, topBinValue)
        topCutIndex=len(smartClipped-1) \
            -int(math.floor(float((100.0-cutLevels[1])/100.0)*len(smartClipped-1)))
        bottomCutIndex=int(math.ceil(float((100.0-cutLevels[1])/100.0)*len(smartClipped-1)))
        topCut=smartClipped[topCutIndex]
        bottomCut=smartClipped[bottomCutIndex]
        anorm=pylab.Normalize(bottomCut, topCut)
    else:
        
        # Normalise using given cut levels
        anorm=pylab.Normalize(cutLevels[0], cutLevels[1])
    
    if cutLevels[0]=="histEq":
        return {'image': imageData.copy()}
    else:
        return {'image': imageData.copy(), 'norm': anorm}

#---------------------------------------------------------------------------------------------------
def resampleToTanProjection(imageData, imageWCS, outputPixDimensions=[600, 600]):
    """Resamples an image and WCS to a tangent plane projection. Purely for plotting purposes
    (e.g., ensuring RA, dec. coordinate axes perpendicular).
    
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS
    @param imageWCS: astWCS.WCS object
    @type outputPixDimensions: list
    @param outputPixDimensions: [width, height] of output image in pixels
    @rtype: dictionary
    @return: image data (np array), updated astWCS WCS object for image, in format {'data', 'wcs'}.
    
    """
    
    RADeg, decDeg=imageWCS.getCentreWCSCoords()
    xPixelScale=imageWCS.getXPixelSizeDeg()
    yPixelScale=imageWCS.getYPixelSizeDeg()
    xSizeDeg, ySizeDeg=imageWCS.getFullSizeSkyDeg()
    xSizePix=int(round(outputPixDimensions[0]))
    ySizePix=int(round(outputPixDimensions[1]))
    xRefPix=xSizePix/2.0
    yRefPix=ySizePix/2.0
    xOutPixScale=xSizeDeg/xSizePix
    yOutPixScale=ySizeDeg/ySizePix
    newHead=pyfits.Header()
    newHead['NAXIS']=2
    newHead['NAXIS1']=xSizePix
    newHead['NAXIS2']=ySizePix
    newHead['CTYPE1']='RA---TAN'
    newHead['CTYPE2']='DEC--TAN'
    newHead['CRVAL1']=RADeg
    newHead['CRVAL2']=decDeg
    newHead['CRPIX1']=xRefPix+1
    newHead['CRPIX2']=yRefPix+1
    newHead['CDELT1']=-xOutPixScale
    newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
    newHead['CUNIT1']='DEG'
    newHead['CUNIT2']='DEG'
    newWCS=astWCS.WCS(newHead, mode='pyfits')
    newImage=np.zeros([ySizePix, xSizePix])

    tanImage=resampleToWCS(newImage, newWCS, imageData, imageWCS, highAccuracy=True, 
                            onlyOverlapping=False)
    
    return tanImage 
    
#---------------------------------------------------------------------------------------------------
def resampleToWCS(im1Data, im1WCS, im2Data, im2WCS, highAccuracy = False, onlyOverlapping = True):
    """Resamples data corresponding to second image (with data im2Data, WCS im2WCS) onto the WCS 
    of the first image (im1Data, im1WCS). The output, resampled image is of the pixel same 
    dimensions of the first image. This routine is for assisting in plotting - performing 
    photometry on the output is not recommended. 
    
    Set highAccuracy == True to sample every corresponding pixel in each image; otherwise only
    every nth pixel (where n is the ratio of the image scales) will be sampled, with values
    in between being set using a linear interpolation (much faster).
    
    Set onlyOverlapping == True to speed up resampling by only resampling the overlapping
    area defined by both image WCSs.
    
    @type im1Data: np array
    @param im1Data: image data array for first image
    @type im1WCS: astWCS.WCS
    @param im1WCS: astWCS.WCS object corresponding to im1Data
    @type im2Data: np array
    @param im2Data: image data array for second image (to be resampled to match first image)
    @type im2WCS: astWCS.WCS
    @param im2WCS: astWCS.WCS object corresponding to im2Data
    @type highAccuracy: bool
    @param highAccuracy: if True, sample every corresponding pixel in each image; otherwise, sample
        every nth pixel, where n = the ratio of the image scales.
    @type onlyOverlapping: bool
    @param onlyOverlapping: if True, only consider the overlapping area defined by both image WCSs
        (speeds things up)
    @rtype: dictionary
    @return: np image data array and associated WCS in format {'data', 'wcs'}
    
    """
    
    resampledData=np.zeros(im1Data.shape)
    
    # Find overlap - speed things up
    # But have a border so as not to require the overlap to be perfect
    # There's also no point in oversampling image 1 if it's much higher res than image 2
    xPixRatio=(im2WCS.getXPixelSizeDeg()/im1WCS.getXPixelSizeDeg())/2.0
    yPixRatio=(im2WCS.getYPixelSizeDeg()/im1WCS.getYPixelSizeDeg())/2.0
    xBorder=xPixRatio*10.0
    yBorder=yPixRatio*10.0
    if highAccuracy == False:
        if xPixRatio > 1:
            xPixStep=int(math.ceil(xPixRatio))
        else:
            xPixStep=1
        if yPixRatio > 1:
            yPixStep=int(math.ceil(yPixRatio))
        else:
            yPixStep=1
    else:
        xPixStep=1
        yPixStep=1
    
    if onlyOverlapping == True:
        overlap=astWCS.findWCSOverlap(im1WCS, im2WCS)
        xOverlap=[overlap['wcs1Pix'][0], overlap['wcs1Pix'][1]]
        yOverlap=[overlap['wcs1Pix'][2], overlap['wcs1Pix'][3]]
        xOverlap.sort()
        yOverlap.sort()
        xMin=int(math.floor(xOverlap[0]-xBorder))
        xMax=int(math.ceil(xOverlap[1]+xBorder))
        yMin=int(math.floor(yOverlap[0]-yBorder))
        yMax=int(math.ceil(yOverlap[1]+yBorder))
        xRemainder=(xMax-xMin) % xPixStep
        yRemainder=(yMax-yMin) % yPixStep
        if xRemainder != 0:
            xMax=xMax+xRemainder
        if yRemainder != 0:
            yMax=yMax+yRemainder
        # Check that we're still within the image boundaries, to be on the safe side
        if xMin < 0:
            xMin=0
        if xMax > im1Data.shape[1]:
            xMax=im1Data.shape[1]
        if yMin < 0:
            yMin=0
        if yMax > im1Data.shape[0]:
            yMax=im1Data.shape[0]
    else:
        xMin=0
        xMax=im1Data.shape[1]
        yMin=0
        yMax=im1Data.shape[0]
    
    for x in range(xMin, xMax, xPixStep):
        for y in range(yMin, yMax, yPixStep):
            RA, dec=im1WCS.pix2wcs(x, y)
            x2, y2=im2WCS.wcs2pix(RA, dec)
            x2=int(round(x2))
            y2=int(round(y2))
            if x2 >= 0 and x2 < im2Data.shape[1] and y2 >= 0 and y2 < im2Data.shape[0]:
                resampledData[y][x]=im2Data[y2][x2]

    # linear interpolation
    if highAccuracy == False:
        for row in range(resampledData.shape[0]):
            vals=resampledData[row, np.arange(xMin, xMax, xPixStep)]
            index2data=interpolate.interp1d(np.arange(0, vals.shape[0], 1), vals)
            interpedVals=index2data(np.arange(0, vals.shape[0]-1, 1.0/xPixStep))
            resampledData[row, xMin:xMin+interpedVals.shape[0]]=interpedVals
        for col in range(resampledData.shape[1]):
            vals=resampledData[np.arange(yMin, yMax, yPixStep), col]
            index2data=interpolate.interp1d(np.arange(0, vals.shape[0], 1), vals)
            interpedVals=index2data(np.arange(0, vals.shape[0]-1, 1.0/yPixStep))
            resampledData[yMin:yMin+interpedVals.shape[0], col]=interpedVals
        
    # Note: should really just copy im1WCS keywords into im2WCS and return that
    # Only a problem if we're using this for anything other than plotting
    return {'data': resampledData, 'wcs': im1WCS.copy()}
    
#---------------------------------------------------------------------------------------------------
def generateContourOverlay(backgroundImageData, backgroundImageWCS, contourImageData, contourImageWCS, \
                            contourLevels, contourSmoothFactor = 0, highAccuracy = False):
    """Rescales an image array to be used as a contour overlay to have the same dimensions as the 
    background image, and generates a set of contour levels. The image array from which the contours 
    are to be generated will be resampled to the same dimensions as the background image data, and 
    can be optionally smoothed using a Gaussian filter. The sigma of the Gaussian filter 
    (contourSmoothFactor) is specified in arcsec.
    
    @type backgroundImageData: np array
    @param backgroundImageData: background image data array
    @type backgroundImageWCS: astWCS.WCS
    @param backgroundImageWCS: astWCS.WCS object of the background image data array
    @type contourImageData: np array
    @param contourImageData: image data array from which contours are to be generated
    @type contourImageWCS: astWCS.WCS
    @param contourImageWCS: astWCS.WCS object corresponding to contourImageData
    @type contourLevels: list
    @param contourLevels: sets the contour levels - available options:
        - values: contourLevels=[list of values specifying each level]
        - linear spacing: contourLevels=['linear', min level value, max level value, number
        of levels] - can use "min", "max" to automatically set min, max levels from image data
        - log spacing: contourLevels=['log', min level value, max level value, number of
        levels] - can use "min", "max" to automatically set min, max levels from image data
    @type contourSmoothFactor: float
    @param contourSmoothFactor: standard deviation (in arcsec) of Gaussian filter for
    pre-smoothing of contour image data (set to 0 for no smoothing)
    @type highAccuracy: bool
    @param highAccuracy: if True, sample every corresponding pixel in each image; otherwise, sample
        every nth pixel, where n = the ratio of the image scales.
    
    """	
    
    # For compromise between speed and accuracy, scale a copy of the background
    # image down to a scale that is one pixel = 1/5 of a pixel in the contour image
    # But only do this if it has CDij keywords as we know how to scale those
    if ("CD1_1" in backgroundImageWCS.header) == True:
        xScaleFactor=backgroundImageWCS.getXPixelSizeDeg()/(contourImageWCS.getXPixelSizeDeg()/5.0)
        yScaleFactor=backgroundImageWCS.getYPixelSizeDeg()/(contourImageWCS.getYPixelSizeDeg()/5.0)
        scaledBackground=scaleImage(backgroundImageData, backgroundImageWCS, (xScaleFactor, yScaleFactor))
        scaled=resampleToWCS(scaledBackground['data'], scaledBackground['wcs'], 
                                contourImageData, contourImageWCS, highAccuracy = highAccuracy)
        scaledContourData=scaled['data']
        scaledContourWCS=scaled['wcs']
        scaledBackground=True
    else:
        scaled=resampleToWCS(backgroundImageData, backgroundImageWCS, 
                                contourImageData, contourImageWCS, highAccuracy = highAccuracy)
        scaledContourData=scaled['data']
        scaledContourWCS=scaled['wcs']
        scaledBackground=False

    if contourSmoothFactor != None and contourSmoothFactor > 0:
        sigmaPix=(contourSmoothFactor/3600.0)/scaledContourWCS.getPixelSizeDeg()
        scaledContourData=ndimage.gaussian_filter(scaledContourData, sigmaPix)
    
    # Various ways of setting the contour levels
    # If just a list is passed in, use those instead
    if contourLevels[0] == "linear":
        if contourLevels[1] == "min":
            xMin=contourImageData.flatten().min()
        else:
            xMin=float(contourLevels[1])
        if contourLevels[2] == "max":
            xMax=contourImageData.flatten().max()
        else:
            xMax=float(contourLevels[2])        
        nLevels=contourLevels[3]
        xStep=(xMax-xMin)/(nLevels-1)
        cLevels=[]
        for j in range(nLevels+1):
            level=xMin+j*xStep
            cLevels.append(level)
    
    elif contourLevels[0] == "log":
        if contourLevels[1] == "min":
            xMin=contourImageData.flatten().min()
        else:
            xMin=float(contourLevels[1])
        if contourLevels[2] == "max":
            xMax=contourImageData.flatten().max()
        else:
            xMax=float(contourLevels[2])     
        if xMin <= 0.0:
            raise Exception("minimum contour level set to <= 0 and log scaling chosen.")
        xLogMin=math.log10(xMin)
        xLogMax=math.log10(xMax)
        nLevels=contourLevels[3]
        xLogStep=(xLogMax-xLogMin)/(nLevels-1)
        cLevels=[]
        prevLevel=0
        for j in range(nLevels+1):
            level=math.pow(10, xLogMin+j*xLogStep)
            cLevels.append(level)			
        
    else:
        cLevels=contourLevels
    
    # Now blow the contour image data back up to the size of the original image   
    if scaledBackground == True:
        scaledBack=scaleImage(scaledContourData, scaledContourWCS, (1.0/xScaleFactor, 1.0/yScaleFactor))['data']
    else:
        scaledBack=scaledContourData
    
    return {'scaledImage': scaledBack, 'contourLevels': cLevels}
    
#---------------------------------------------------------------------------------------------------
def saveBitmap(outputFileName, imageData, cutLevels, size, colorMapName):
    """Makes a bitmap image from an image array; the image format is specified by the
    filename extension. (e.g. ".jpg" =JPEG, ".png"=PNG).
    
    @type outputFileName: string
    @param outputFileName: filename of output bitmap image
    @type imageData: np array
    @param imageData: image data array
    @type cutLevels: list
    @param cutLevels: sets the image scaling - available options:
        - pixel values: cutLevels=[low value, high value].
        - histogram equalisation: cutLevels=["histEq", number of bins ( e.g. 1024)]
        - relative: cutLevels=["relative", cut per cent level (e.g. 99.5)]
        - smart: cutLevels=["smart", cut per cent level (e.g. 99.5)]
    ["smart", 99.5] seems to provide good scaling over a range of different images. 
    @type size: int
    @param size: size of output image in pixels
    @type colorMapName: string
    @param colorMapName: name of a standard matplotlib colormap, e.g. "hot", "cool", "gray"
    etc. (do "help(pylab.colormaps)" in the Python interpreter to see available options)
    
    """		
    
    cut=intensityCutImage(imageData, cutLevels)
    
    # Make plot
    aspectR=float(cut['image'].shape[0])/float(cut['image'].shape[1])
    pylab.figure(figsize=(10,10*aspectR))
    pylab.axes([0,0,1,1])
        
    try:
        colorMap=pylab.cm.get_cmap(colorMapName)
    except AssertionError:
        raise Exception(colorMapName+" is not a defined matplotlib colormap.")
    
    if cutLevels[0]=="histEq":
        pylab.imshow(cut['image'],  interpolation="bilinear", origin='lower', cmap=colorMap)
    
    else:
        pylab.imshow(cut['image'],  interpolation="bilinear",  norm=cut['norm'], origin='lower',
            cmap=colorMap)

    pylab.axis("off")
    
    pylab.savefig("out_astImages.png")	
    pylab.close("all")
    
    try:
        from PIL import Image
    except:
        raise Exception("astImages.saveBitmap requires the Python Imaging Library to be installed.")
    im=Image.open("out_astImages.png")
    im.thumbnail((int(size),int(size)))
    im.save(outputFileName)
    
    os.remove("out_astImages.png")

#---------------------------------------------------------------------------------------------------
def saveContourOverlayBitmap(outputFileName, backgroundImageData, backgroundImageWCS, cutLevels, \
                                size, colorMapName, contourImageData, contourImageWCS, \
                                contourSmoothFactor, contourLevels, contourColor, contourWidth):
    """Makes a bitmap image from an image array, with a set of contours generated from a
    second image array overlaid. The image format is specified by the file extension
    (e.g. ".jpg"=JPEG, ".png"=PNG). The image array from which the contours are to be generated
    can optionally be pre-smoothed using a Gaussian filter. 
    
    @type outputFileName: string
    @param outputFileName: filename of output bitmap image
    @type backgroundImageData: np array
    @param backgroundImageData: background image data array
    @type backgroundImageWCS: astWCS.WCS
    @param backgroundImageWCS: astWCS.WCS object of the background image data array
    @type cutLevels: list
    @param cutLevels: sets the image scaling - available options:
        - pixel values: cutLevels=[low value, high value].
        - histogram equalisation: cutLevels=["histEq", number of bins ( e.g. 1024)]
        - relative: cutLevels=["relative", cut per cent level (e.g. 99.5)]
        - smart: cutLevels=["smart", cut per cent level (e.g. 99.5)]
    ["smart", 99.5] seems to provide good scaling over a range of different images. 
    @type size: int
    @param size: size of output image in pixels
    @type colorMapName: string
    @param colorMapName: name of a standard matplotlib colormap, e.g. "hot", "cool", "gray"
    etc. (do "help(pylab.colormaps)" in the Python interpreter to see available options)
    @type contourImageData: np array
    @param contourImageData: image data array from which contours are to be generated
    @type contourImageWCS: astWCS.WCS
    @param contourImageWCS: astWCS.WCS object corresponding to contourImageData
    @type contourSmoothFactor: float
    @param contourSmoothFactor: standard deviation (in pixels) of Gaussian filter for
    pre-smoothing of contour image data (set to 0 for no smoothing)
    @type contourLevels: list
    @param contourLevels: sets the contour levels - available options:
        - values: contourLevels=[list of values specifying each level]
        - linear spacing: contourLevels=['linear', min level value, max level value, number
        of levels] - can use "min", "max" to automatically set min, max levels from image data
        - log spacing: contourLevels=['log', min level value, max level value, number of
        levels] - can use "min", "max" to automatically set min, max levels from image data
    @type contourColor: string
    @param contourColor: color of the overlaid contours, specified by the name of a standard
    matplotlib color, e.g., "black", "white", "cyan"
    etc. (do "help(pylab.colors)" in the Python interpreter to see available options)
    @type contourWidth: int
    @param contourWidth: width of the overlaid contours
    
    """	
    
    cut=intensityCutImage(backgroundImageData, cutLevels)
    
    # Make plot of just the background image
    aspectR=float(cut['image'].shape[0])/float(cut['image'].shape[1])
    pylab.figure(figsize=(10,10*aspectR))
    pylab.axes([0,0,1,1])
        
    try:
        colorMap=pylab.cm.get_cmap(colorMapName)
    except AssertionError:
        raise Exception(colorMapName+" is not a defined matplotlib colormap.")
    
    if cutLevels[0]=="histEq":
        pylab.imshow(cut['image'],  interpolation="bilinear", origin='lower', cmap=colorMap)
    
    else:
        pylab.imshow(cut['image'],  interpolation="bilinear",  norm=cut['norm'], origin='lower',
            cmap=colorMap)

    pylab.axis("off")

    # Add the contours
    contourData=generateContourOverlay(backgroundImageData, backgroundImageWCS, contourImageData, \
                                        contourImageWCS, contourLevels, contourSmoothFactor)
    
    pylab.contour(contourData['scaledImage'], contourData['contourLevels'], colors=contourColor,
        linewidths=contourWidth)	
            
    pylab.savefig("out_astImages.png")	
    pylab.close("all")
    
    try:
        from PIL import Image
    except:
        raise Exception("astImages.saveContourOverlayBitmap requires the Python Imaging Library to be installed")
    
    im=Image.open("out_astImages.png")
    im.thumbnail((int(size),int(size)))
    im.save(outputFileName)
        
    os.remove("out_astImages.png")
    
#---------------------------------------------------------------------------------------------------
def saveFITS(outputFileName, imageData, imageWCS = None):
    """Writes an image array to a new .fits file.
    
    @type outputFileName: string
    @param outputFileName: filename of output FITS image
    @type imageData: np array
    @param imageData: image data array
    @type imageWCS: astWCS.WCS object
    @param imageWCS: image WCS object
    
    @note: If imageWCS=None, the FITS image will be written with a rudimentary header containing
    no meta data.
    
    """
    
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    
    # There a fudge here for handling both pyfits and astropy.io.fits headers
    # Removed from version 0.10.0+ (supporting astropy only)
    if imageWCS != None:
        hdu=pyfits.PrimaryHDU(None, imageWCS.header)
    else:
        hdu=pyfits.PrimaryHDU(None, None)
    
    newImg=pyfits.HDUList()
    hdu.data=imageData
    newImg.append(hdu)
    newImg.writeto(outputFileName)
    newImg.close()
    
#---------------------------------------------------------------------------------------------------
def histEq(inputArray, numBins):
    """Performs histogram equalisation of the input np array.
    
    @type inputArray: np array
    @param inputArray: image data array
    @type numBins: int
    @param numBins: number of bins in which to perform the operation (e.g. 1024)
    @rtype: np array
    @return: image data array
    
    """
    
    imageData=inputArray
    
    # histogram equalisation: we want an equal number of pixels in each intensity range
    sortedDataIntensities=np.sort(np.ravel(imageData))	
    median=np.median(sortedDataIntensities)
    
    # Make cumulative histogram of data values, simple min-max used to set bin sizes and range
    dataCumHist=np.zeros(numBins)
    minIntensity=sortedDataIntensities.min()	
    maxIntensity=sortedDataIntensities.max()
    histRange=maxIntensity-minIntensity
    binWidth=histRange/float(numBins-1)
    for i in range(len(sortedDataIntensities)):
        binNumber=int(math.ceil((sortedDataIntensities[i]-minIntensity)/binWidth))
        addArray=np.zeros(numBins)
        onesArray=np.ones(numBins-binNumber)
        onesRange=list(range(binNumber, numBins))
        np.put(addArray, onesRange, onesArray)
        dataCumHist=dataCumHist+addArray
                
    # Make ideal cumulative histogram
    idealValue=dataCumHist.max()/float(numBins)
    idealCumHist=np.arange(idealValue, dataCumHist.max()+idealValue, idealValue)
    
    # Map the data to the ideal
    for y in range(imageData.shape[0]):
        for x in range(imageData.shape[1]):
            # Get index corresponding to dataIntensity
            intensityBin=int(math.ceil((imageData[y][x]-minIntensity)/binWidth))
            
            # Guard against rounding errors (happens rarely I think)
            if intensityBin<0:
                intensityBin=0
            if intensityBin>len(dataCumHist)-1:
                intensityBin=len(dataCumHist)-1
        
            # Get the cumulative frequency corresponding intensity level in the data
            dataCumFreq=dataCumHist[intensityBin]
            
            # Get the index of the corresponding ideal cumulative frequency
            idealBin=np.searchsorted(idealCumHist, dataCumFreq)
            idealIntensity=(idealBin*binWidth)+minIntensity
            imageData[y][x]=idealIntensity	
        
    return imageData

#---------------------------------------------------------------------------------------------------
def normalise(inputArray, clipMinMax):
    """Clips the inputArray in intensity and normalises the array such that minimum and maximum
    values are 0, 1. Clip in intensity is specified by clipMinMax, a list in the format 
    [clipMin, clipMax]
    
    Used for normalising image arrays so that they can be turned into RGB arrays that matplotlib
    can plot (see L{astPlots.ImagePlot}).
    
    @type inputArray: np array
    @param inputArray: image data array
    @type clipMinMax: list
    @param clipMinMax: [minimum value of clipped array, maximum value of clipped array]
    @rtype: np array
    @return: normalised array with minimum value 0, maximum value 1

    """
    clipped=inputArray.clip(clipMinMax[0], clipMinMax[1])
    slope=1.0/(clipMinMax[1]-clipMinMax[0])
    intercept=-clipMinMax[0]*slope
    clipped=clipped*slope+intercept
    
    return clipped
    
