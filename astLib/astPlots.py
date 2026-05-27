"""Module for producing astronomical plots.

(c) 2007-2024 Matt Hilton

This module provides the matplotlib powered ImagePlot class, which is designed to be flexible.
ImagePlots can have RA, Dec. coordinate axes, contour overlays, and have objects marked in them,
using WCS coordinates. RGB plots are supported too.

"""

import math
from . import astImages
from . import astWCS
from . import astCoords
import numpy as np
import astropy.io.fits as pyfits
from scipy import interpolate
import pylab
import matplotlib.patches as patches
import sys

# Handle unicode python 2 and 3
if sys.version < '3':
    import codecs
    def u(x):
        return codecs.unicode_escape_decode(x)[0]
else:
    def u(x):
        return x

DEC_TICK_STEPS=[{'deg': 1.0/60.0/60.0,  'unit': "s"},
                {'deg': 2.0/60.0/60.0,  'unit': "s"},
                {'deg': 5.0/60.0/60.0,  'unit': "s"},
                {'deg': 10.0/60.0/60.0, 'unit': "s"},
                {'deg': 30.0/60.0/60.0, 'unit': "s"},
                {'deg': 1.0/60.0,       'unit': "m"},
                {'deg': 2.0/60.0,       'unit': "m"},
                {'deg': 5.0/60.0,       'unit': "m"},
                {'deg': 15.0/60.0,      'unit': "m"},
                {'deg': 30.0/60.0,      'unit': "m"},
                {'deg': 1.0,            'unit': "d"},
                {'deg': 2.0,            'unit': "d"},
                {'deg': 4.0,            'unit': "d"},
                {'deg': 5.0,            'unit': "d"},
                {'deg': 10.0,           'unit': "d"},
                {'deg': 20.0,           'unit': "d"},
                {'deg': 30.0,           'unit': "d"}]
"""Defines the possible coordinate label steps on the delination axis in
sexagesimal mode. Dictionary format: {'deg', 'unit'}"""

RA_TICK_STEPS=[ {'deg': (0.5/60.0/60.0/24.0)*360.0,  'unit': "s"},
                {'deg': (1.0/60.0/60.0/24.0)*360.0,  'unit': "s"},
                {'deg': (2.0/60.0/60.0/24.0)*360.0,  'unit': "s"},
                {'deg': (4.0/60.0/60.0/24.0)*360.0,  'unit': "s"},
                {'deg': (5.0/60.0/60.0/24.0)*360.0,  'unit': "s"},
                {'deg': (10.0/60.0/60.0/24.0)*360.0, 'unit': "s"},
                {'deg': (20.0/60.0/60.0/24.0)*360.0, 'unit': "s"},
                {'deg': (30.0/60.0/60.0/24.0)*360.0, 'unit': "s"},
                {'deg': (1.0/60.0/24.0)*360.0,       'unit': "m"},
                {'deg': (2.0/60.0/24.0)*360.0,       'unit': "m"},
                {'deg': (5.0/60.0/24.0)*360.0,       'unit': "m"},
                {'deg': (10.0/60.0/24.0)*360.0,      'unit': "m"},
                {'deg': (20.0/60.0/24.0)*360.0,      'unit': "m"},
                {'deg': (30.0/60.0/24.0)*360.0,      'unit': "m"},
                {'deg': (1.0/24.0)*360.0,            'unit': "h"},
                {'deg': (3.0/24.0)*360.0,            'unit': "h"},
                {'deg': (6.0/24.0)*360.0,            'unit': "h"},
                {'deg': (12.0/24.0)*360.0,           'unit': "h"}]
"""Defines the possible coordinate label steps on the right ascension axis in
sexagesimal mode. Dictionary format: {'deg', 'unit'}"""

DECIMAL_TICK_STEPS=[0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 30.0, 90.0]
"""Defines the possible coordinate label steps on both coordinate axes in decimal degrees mode."""

DEG = u("\N{DEGREE SIGN}")

PRIME = "$^\prime$"

DOUBLE_PRIME = "$^{\prime\prime}$"

#---------------------------------------------------------------------------------------------------
class ImagePlot:
    """This class describes a matplotlib image plot containing an astronomical image with an
    associated WCS.

    Objects within the image boundaries can be marked by passing their WCS coordinates to
    :meth:`~astLib.astPlots.ImagePlot.addPlotObjects`.

    Other images can be overlaid using :meth:`~astLib.astPlots.ImagePlot.addContourOverlay`.

    For images rotated with North at the top, East at the left (as can be done using
    :func:`~astLib.astImages.clipRotatedImageSectionWCS` or :func:`~astLib.astImages.resampleToTanProjection`), WCS coordinate
    axes can be plotted, with tick marks set appropriately for the image size. Otherwise, a compass
    can be plotted showing the directions of North and East in the image.

    RGB images are also supported.

    The plot can of course be tweaked further after creation using matplotlib/pylab commands.

    """
    def __init__(self, imageData, imageWCS, axes = [0.1,0.1,0.8,0.8], \
        cutLevels = ["smart", 99.5], colorMapName = "gray", title = None, axesLabels = "sexagesimal", \
        axesFontFamily="serif", axesFontSize=12.0, RATickSteps="auto", decTickSteps="auto",
        colorBar = False, interpolation = "bilinear"):
        """Makes an ImagePlot from the given image array and astWCS. For coordinate axes to work, the
        image and WCS should have been rotated such that East is at the left, North is at the top
        (see e.g. :func:`~astLib.astImages.clipRotatedImageSectionWCS`, or :func:`~astLib.astImages.resampleToTanProjection`).

        If imageData is given as a list in the format [r, g, b], a color RGB plot will be made. However,
        in this case the cutLevels must be specified manually for each component as a list -
        i.e. cutLevels = [[r min, r max], [g min, g max], [b min, b max]]. In this case of course, the
        colorMap will be ignored. All r, g, b image arrays must have the same dimensions.

        Set axesLabels = None to make a plot without coordinate axes plotted.

        The axes can be marked in either sexagesimal or decimal celestial coordinates. If RATickSteps
        or decTickSteps are set to "auto", the appropriate axis scales will be determined automatically
        from the size of the image array and associated WCS. The tick step sizes can be overidden.
        If the coordinate axes are in sexagesimal format a dictionary in the format {'deg', 'unit'} is
        needed (see :data:`~astLib.astPlots.RA_TICK_STEPS` and :data:`~astLib.astPlots.DEC_TICK_STEPS` for examples). If the coordinate axes are in
        decimal format, the tick step size is specified simply in RA, dec decimal degrees.

        Args:
            imageData (numpy.ndarray or list): image data array or list of np arrays [r, g, b]
            imageWCS (astWCS.WCS): astWCS.WCS object
            axes (list): specifies where in the current figure to draw the finder chart (see pylab.axes)
            cutLevels (list): sets the image scaling - available options:

                - pixel values: ``cutLevels=[low value, high value]``
                - histogram equalisation: ``cutLevels=["histEq", number of bins (e.g. 1024)]``
                - relative: ``cutLevels=["relative", cut per cent level (e.g. 99.5)]``
                - smart: ``cutLevels=["smart", cut per cent level (e.g. 99.5)]``

                ``["smart", 99.5]`` seems to provide good scaling over a range of different images.
                For RGB images, cut levels must be specified manually as a list:
                ``[[r min, r max], [g min, g max], [b min, b max]]``
            colorMapName (str): name of a standard matplotlib colormap, e.g. "hot", "cool", "gray"
            title (str): optional title for the plot
            axesLabels (str): either "sexagesimal" (for H:M:S, D:M:S), "decimal" (for decimal degrees)
                or None (for no coordinate axes labels)
            axesFontFamily (str): matplotlib fontfamily, e.g. 'serif', 'sans-serif' etc.
            axesFontSize (float): font size of axes labels and titles (in points)
            colorBar (bool): if True, plot a vertical color bar at the side of the image indicating
                the intensity scale
            interpolation (str): interpolation to apply to the image plot (see the documentation for
                the matplotlib.pylab.imshow command)

        """

        self.RADeg, self.decDeg=imageWCS.getCentreWCSCoords()
        self.wcs=imageWCS

        # Handle case where imageData is [r, g, b]
        if type(imageData) == list:
            if len(imageData) == 3:
                if len(cutLevels) == 3:
                    r=astImages.normalise(imageData[0], cutLevels[0])
                    g=astImages.normalise(imageData[1], cutLevels[1])
                    b=astImages.normalise(imageData[2], cutLevels[2])
                    rgb=np.array([r.transpose(), g.transpose(), b.transpose()])
                    rgb=rgb.transpose()
                    self.data=rgb
                    self.rgbImage=True
                else:
                    raise Exception("tried to create a RGB array, but cutLevels is not a list of 3 lists")

            else:
                raise Exception("tried to create a RGB array but imageData is not a list of 3 arrays")
        else:
            self.data=imageData
            self.rgbImage=False

        self.axes=pylab.axes(axes)
        self.cutLevels=cutLevels
        self.colorMapName=colorMapName
        self.title=title
        self.axesLabels=axesLabels
        self.colorBar=colorBar
        self.axesFontSize=axesFontSize
        self.axesFontFamily=axesFontFamily

        self.flipXAxis=False
        self.flipYAxis=False

        self.interpolation=interpolation

        if self.axesLabels != None:

            # Allow user to override the automatic coord tick spacing
            if self.axesLabels == "sexagesimal":
                if RATickSteps != "auto":
                    if type(RATickSteps) != dict or "deg" not in list(RATickSteps.keys()) \
                            or "unit" not in list(RATickSteps.keys()):
                        raise Exception("RATickSteps needs to be in format {'deg', 'unit'} for sexagesimal axes labels")
                if decTickSteps != "auto":
                    if type(decTickSteps) != dict or "deg" not in list(decTickSteps.keys()) \
                            or "unit" not in list(decTickSteps.keys()):
                        raise Exception("decTickSteps needs to be in format {'deg', 'unit'} for sexagesimal axes labels")
            elif self.axesLabels == "decimal":
                if RATickSteps != "auto":
                    if type(RATickSteps) != float:
                        raise Exception("RATickSteps needs to be a float (if not 'auto') for decimal axes labels")
                if decTickSteps != "auto":
                    if type(decTickSteps) != float:
                        raise Exception("decTickSteps needs to be a float (if not 'auto') for decimal axes labels")
            self.RATickSteps=RATickSteps
            self.decTickSteps=decTickSteps

            self.calcWCSAxisLabels(axesLabels = self.axesLabels)

        # this list stores objects to overplot, add to it using addPlotObjects()
        self.plotObjects=[]

        # this list stores image data to overlay as contours, add to it using addContourOverlay()
        self.contourOverlays=[]

        self.draw()


    def draw(self):
        """Redraws the ImagePlot.

        """

        pylab.axes(self.axes)
        pylab.cla()

        if self.title != None:
            pylab.title(self.title)
        try:
            colorMap=pylab.cm.get_cmap(self.colorMapName)
        except AssertionError:
            raise Exception(self.colorMapName+"is not a defined matplotlib colormap.")

        if self.rgbImage == False:
            self.cutImage=astImages.intensityCutImage(self.data, self.cutLevels)
            if self.cutLevels[0]=="histEq":
                pylab.imshow(self.cutImage['image'],  interpolation=self.interpolation, origin='lower', cmap=colorMap)
            else:
                pylab.imshow(self.cutImage['image'],  interpolation=self.interpolation,  norm=self.cutImage['norm'], \
                            origin='lower', cmap=colorMap)
        else:
            pylab.imshow(self.data, interpolation="bilinear", origin='lower')

        if self.colorBar == True:
            pylab.colorbar(shrink=0.8)

        for c in self.contourOverlays:
            pylab.contour(c['contourData']['scaledImage'], c['contourData']['contourLevels'],
                            colors=c['color'], linewidths=c['width'])

        for p in self.plotObjects:
            for x, y, l in zip(p['x'], p['y'], p['objLabels']):
                if p['symbol'] == "circle":
                    c=patches.Circle((x, y), radius=p['sizePix']/2.0, fill=p['fill'], color=p['color'],
                                        linewidth=p['width'])
                    self.axes.add_patch(c)
                elif p['symbol'] == "box":
                    c=patches.Rectangle((x-p['sizePix']/2, y-p['sizePix']/2), p['sizePix'], p['sizePix'],
                        fill=p['fill'], color=p['color'], linewidth=p['width'])
                    self.axes.add_patch(c)
                elif p['symbol'] == "cross":
                    pylab.plot([x-p['sizePix']/2, x+p['sizePix']/2], [y, y], linestyle='-',
                        linewidth=p['width'], color= p['color'])
                    pylab.plot([x, x], [y-p['sizePix']/2, y+p['sizePix']/2], linestyle='-',
                        linewidth=p['width'], color= p['color'])
                elif p['symbol'] == "diamond":
                    c=patches.RegularPolygon([x, y], 4, radius=p['sizePix']/2, orientation=0,
                                             color=p['color'], fill=p['fill'], linewidth=p['width'])
                    self.axes.add_patch(c)
                if l != None:
                    pylab.text(x, y+p['sizePix']/1.5, l, horizontalalignment='center', \
                                fontsize=p['objLabelSize'], color=p['color'])

            if p['symbol'] == "compass":
                x=p['x'][0]
                y=p['y'][0]
                ra=p['RA'][0]
                dec=p['dec'][0]

                westPoint,eastPoint,southPoint,northPoint=astCoords.calcRADecSearchBox(ra, dec, p['sizeArcSec']/3600.0/2.0)
                northPix=self.wcs.wcs2pix(ra, northPoint)
                eastPix=self.wcs.wcs2pix(eastPoint, dec)

                edx=eastPix[0]-x
                edy=eastPix[1]-y
                ndx=northPix[0]-x
                ndy=northPix[1]-y
                nArrow=patches.Arrow(x, y, ndx, ndy, edgecolor=p['color'], facecolor=p['color'], width=p['width'])
                eArrow=patches.Arrow(x, y, edx, edy, edgecolor=p['color'], facecolor=p['color'], width=p['width'])
                self.axes.add_patch(nArrow)
                self.axes.add_patch(eArrow)
                pylab.text(x+ndx+ndx*0.2, y+ndy+ndy*0.2, "N", horizontalalignment='center',
                                verticalalignment='center', fontsize=p['objLabelSize'], color=p['color'])
                pylab.text(x+edx+edx*0.2, y+edy+edy*0.2, "E", horizontalalignment='center',
                                verticalalignment='center', fontsize=p['objLabelSize'], color=p['color'])

            if p['symbol'] == "scaleBar":
                x=p['x'][0]
                y=p['y'][0]
                ra=p['RA'][0]
                dec=p['dec'][0]

                westPoint,eastPoint,southPoint,northPoint=astCoords.calcRADecSearchBox(ra, dec, p['sizeArcSec']/3600.0/2.0)
                northPix=self.wcs.wcs2pix(ra, northPoint)
                eastPix=self.wcs.wcs2pix(eastPoint, dec)
                edx=eastPix[0]-x
                edy=eastPix[1]-y
                ndx=northPix[0]-x
                ndy=northPix[1]-y

                if p['style'] == "arrows":
                    eArrow=patches.Arrow(x, y, edx, edy, edgecolor=p['color'], facecolor=p['color'], width=p['width'])
                    wArrow=patches.Arrow(x, y, -edx, edy, edgecolor=p['color'], facecolor=p['color'], width=p['width'])
                    self.axes.add_patch(eArrow)
                    self.axes.add_patch(wArrow)
                elif p['style'] == "whiskers":
                    ewArrow=patches.FancyArrowPatch(posA = (x+edx, y), posB = (x-edx,y+edy), edgecolor=p['color'], facecolor=p['color'], linewidth = p['width'], arrowstyle = '|-|')
                    self.axes.add_patch(ewArrow)

                # Work out label
                if p['scaleBarLabel'] == None:
                    scaleLabel=None
                    if p['sizeArcSec'] < 60.0:
                        scaleLabel="%.0f %s" % (p['sizeArcSec'], DOUBLE_PRIME)
                    elif p['sizeArcSec'] >= 60.0 and p['sizeArcSec'] <  3600.0:
                        scaleLabel="%.0f %s" % (p['sizeArcSec']/60.0, PRIME)
                    else:
                        scaleLabel="%.0f %s" % (p['sizeArcSec']/3600.0, DEG)
                else:
                    scaleLabel=p['scaleBarLabel']
                pylab.text(x, y+0.025*self.data.shape[1], scaleLabel, horizontalalignment='center',
                           verticalalignment='center', fontsize=p['objLabelSize'], color=p['color'])

        if self.axesLabels != None:
            pylab.xticks(self.ticsRA[0], self.ticsRA[1], weight='normal', family=self.axesFontFamily, \
                                    fontsize=self.axesFontSize)
            pylab.yticks(self.ticsDec[0], self.ticsDec[1], weight='normal', family=self.axesFontFamily, \
                                    fontsize=self.axesFontSize)
            pylab.xlabel(self.RAAxisLabel, family=self.axesFontFamily, fontsize=self.axesFontSize)
            pylab.ylabel(self.decAxisLabel, family=self.axesFontFamily, fontsize=self.axesFontSize)
        else:
            pylab.xticks([], [])
            pylab.yticks([], [])
            pylab.xlabel("")
            pylab.ylabel("")

        if self.flipXAxis == False:
            pylab.xlim(0, self.data.shape[1]-1)
        else:
            pylab.xlim(self.data.shape[1]-1, 0)
        if self.flipYAxis == False:
            pylab.ylim(0, self.data.shape[0]-1)
        else:
            pylab.ylim(self.data.shape[0]-1, 0)


    def addContourOverlay(self, contourImageData, contourWCS, tag, levels = ["linear", "min", "max", 5],
                             width = 1, color = "white", smooth = 0, highAccuracy = False):
        """Adds image data to the ImagePlot as a contour overlay. The contours can be removed using
        :meth:`~astLib.astPlots.ImagePlot.removeContourOverlay`. If a contour overlay already exists with this tag, it will be replaced.

        Args:
            contourImageData (numpy.ndarray): image data array from which contours are to be generated
            contourWCS (astWCS.WCS): astWCS.WCS object for the image to be contoured
            tag (str): identifying tag for this set of contours
            levels (list): sets the contour levels - available options:

                - values: ``contourLevels=[list of values specifying each level]``
                - linear spacing: ``contourLevels=['linear', min level value, max level value, number of levels]``
                  (use "min", "max" to automatically set min, max levels from image data)
                - log spacing: ``contourLevels=['log', min level value, max level value, number of levels]``
                  (use "min", "max" to automatically set min, max levels from image data)

            width (int): width of the overlaid contours
            color (str): color of the overlaid contours, specified by the name of a standard matplotlib
                color, e.g., "black", "white", "cyan"
            smooth (float): standard deviation (in arcsec) of Gaussian filter for pre-smoothing of
                contour image data (set to 0 for no smoothing)
            highAccuracy (bool): if True, sample every corresponding pixel in each image; otherwise,
                sample every nth pixel, where n = the ratio of the image scales

        """

        if self.rgbImage == True:
            backgroundData=self.data[:,:,0]
        else:
            backgroundData=self.data
        contourData=astImages.generateContourOverlay(backgroundData, self.wcs, contourImageData, \
                                              contourWCS, levels, smooth, highAccuracy = highAccuracy)

        alreadyGot=False
        for c in self.contourOverlays:
            if c['tag'] == tag:
                c['contourData']=contourData
                c['tag']=tag
                c['color']=color
                c['width']=width
                alreadyGot=True

        if alreadyGot == False:
            self.contourOverlays.append({'contourData': contourData, 'tag': tag, 'color': color, \
                                            'width': width})
        self.draw()


    def removeContourOverlay(self, tag):
        """Removes the contourOverlay from the ImagePlot corresponding to the tag.

        Args:
            tag (str): tag for contour overlay in ImagePlot.contourOverlays to be removed

        """

        index=0
        for p in self.contourOverlays:
            if p['tag'] == tag:
                self.plotObjects.remove(self.plotObjects[index])
            index=index+1
        self.draw()


    def addPlotObjects(self, objRAs, objDecs, tag, symbol="circle", size=4.0, width=1.0, color="yellow",\
                       fill = False, objLabels = None, objLabelSize = 12.0):
        """Add objects with RA, dec coords objRAs, objDecs to the ImagePlot. Only objects that fall within
        the image boundaries will be plotted.

        symbol specifies the type of symbol with which to mark the object in the image. The following
        values are allowed:

        - "circle"
        - "box"
        - "cross"
        - "diamond"

        size specifies the diameter in arcsec of the symbol (if plotSymbol == "circle"), or the width
        of the box in arcsec (if plotSymbol == "box")

        width specifies the thickness of the symbol lines in pixels

        color can be any valid matplotlib color (e.g. "red", "green", etc.)

        The objects can be removed from the plot by using :meth:`~astLib.astPlots.ImagePlot.removePlotObjects`, and then calling
        draw(). If the ImagePlot already has a set of plotObjects with the same tag, they will be
        replaced.

        Args:
            objRAs (numpy.ndarray or list): object RA coords in decimal degrees
            objDecs (numpy.ndarray or list): corresponding object Dec. coords in decimal degrees
            tag (str): identifying tag for this set of objects
            symbol (str): either "circle", "box", "cross", or "diamond"
            size (float): size of symbols to plot (radius in arcsec, or width of box)
            width (float): width of symbols in pixels
            color (str): any valid matplotlib color string, e.g. "red", "green" etc.
            fill (bool): if True, fill symbols
            objLabels (list): text labels to plot next to objects in figure
            objLabelSize (float): size of font used for object labels (in points)

        """

        pixCoords=self.wcs.wcs2pix(objRAs, objDecs)

        xMax=self.data.shape[1]
        yMax=self.data.shape[0]

        if objLabels is None:
            objLabels=[None]*len(objRAs)

        xInPlot=[]
        yInPlot=[]
        RAInPlot=[]
        decInPlot=[]
        labelInPlot=[]
        for p, r, d, l in zip(pixCoords, objRAs, objDecs, objLabels):
            if p[0] >= 0 and p[0] < xMax and p[1] >= 0 and p[1] < yMax:
                xInPlot.append(p[0])
                yInPlot.append(p[1])
                RAInPlot.append(r)
                decInPlot.append(d)
                labelInPlot.append(l)

        xInPlot=np.array(xInPlot)
        yInPlot=np.array(yInPlot)
        RAInPlot=np.array(RAInPlot)
        decInPlot=np.array(decInPlot)

        # Size of symbols in pixels in plot - converted from arcsec
        sizePix=(size/3600.0)/self.wcs.getPixelSizeDeg()

        alreadyGot=False
        for p in self.plotObjects:
            if p['tag'] == tag:
                p['x']=xInPlot
                p['y']=yInPlot
                p['RA']=RAInPlot
                p['dec']=decInPlot
                p['tag']=tag
                p['objLabels']=objLabels
                p['symbol']=symbol
                p['sizePix']=sizePix
                p['sizeArcSec']=size
                p['width']=width
                p['color']=color
                p['objLabelSize']=objLabelSize
                p['fill']=fill
                alreadyGot=True

        if alreadyGot == False:
            self.plotObjects.append({'x': xInPlot, 'y': yInPlot, 'RA': RAInPlot, 'dec': decInPlot,
                                'tag': tag, 'objLabels': labelInPlot, 'symbol': symbol,
                                'sizePix': sizePix, 'width': width, 'color': color,
                                'objLabelSize': objLabelSize, 'sizeArcSec': size, 'fill': fill})
        self.draw()


    def removePlotObjects(self, tag):
        """Removes the plotObjects from the ImagePlot corresponding to the tag. The plot must be redrawn
        for the change to take effect.

        Args:
            tag (str): tag for set of objects in ImagePlot.plotObjects to be removed

        """

        index=0
        for p in self.plotObjects:
            if p['tag'] == tag:
                self.plotObjects.remove(self.plotObjects[index])
            index=index+1
        self.draw()


    def addCompass(self, location, sizeArcSec, color = "white", fontSize = 12, \
                        width = 20.0):
        """Adds a compass to the ImagePlot at the given location ('N', 'NE', 'E', 'SE', 'S',
        'SW', 'W', or 'NW'). Note these aren't directions on the WCS coordinate grid, they are
        relative positions on the plot - so N is top centre, NE is top right, SW is bottom right etc..
        Alternatively, pixel coordinates (x, y) in the image can be given.

        Args:
            location (str or tuple): location in the plot where the compass is drawn:

                - string: N, NE, E, SE, S, SW, W or NW
                - tuple: (x, y)

            sizeArcSec (float): length of the compass arrows on the plot in arc seconds
            color (str): any valid matplotlib color string
            fontSize (float): size of font used to label N and E, in points
            width (float): width of arrows used to mark compass

        """

        if type(location) == str:
            cRADeg, cDecDeg=self.wcs.getCentreWCSCoords()
            RAMin, RAMax, decMin, decMax=self.wcs.getImageMinMaxWCSCoords()
            westPoint,eastPoint,southPoint,northPoint=astCoords.calcRADecSearchBox(cRADeg, cDecDeg, sizeArcSec/3600.0/2.0)
            sizeRADeg=eastPoint-westPoint
            sizeDecDeg=northPoint-southPoint
            xSizePix=(sizeArcSec/3600.0)/self.wcs.getXPixelSizeDeg()
            ySizePix=(sizeArcSec/3600.0)/self.wcs.getYPixelSizeDeg()
            X=self.data.shape[1]
            Y=self.data.shape[0]
            xBufferPix=0.5*xSizePix
            yBufferPix=0.5*ySizePix
            cx, cy=self.wcs.wcs2pix(cRADeg, cDecDeg)
            foundLocation=False
            x=cy
            y=cx
            if self.wcs.isFlipped() == False:
                if location.find("N") != -1:
                    y=Y-2*yBufferPix
                    foundLocation=True
                if location.find("S") != -1:
                    y=yBufferPix
                    foundLocation=True
                if location.find("E") != -1:
                    x=xBufferPix*2
                    foundLocation=True
                if location.find("W") != -1:
                    x=X-xBufferPix
                    foundLocation=True
            else:
                if location.find("S") != -1:
                    y=Y-2*yBufferPix
                    foundLocation=True
                if location.find("N") != -1:
                    y=yBufferPix
                    foundLocation=True
                if location.find("W") != -1:
                    x=xBufferPix*2
                    foundLocation=True
                if location.find("E") != -1:
                    x=X-xBufferPix
                    foundLocation=True
            if foundLocation == False:
                raise Exception("didn't understand location string for scale bar (should be e.g. N, S, E, W).")
            RADeg, decDeg=self.wcs.pix2wcs(x, y)
        elif type(location) == tuple or type(location) == list:
            x, y=location
            RADeg, decDeg=self.wcs.pix2wcs(x, y)
        else:
            raise Exception("didn't understand location for scale bar - should be string or tuple.")

        alreadyGot=False
        for p in self.plotObjects:
            if p['tag'] == "compass":
                p['x']=[x]
                p['y']=[y]
                p['RA']=[RADeg]
                p['dec']=[decDeg]
                p['tag']="compass"
                p['objLabels']=[None]
                p['symbol']="compass"
                p['sizeArcSec']=sizeArcSec
                p['width']=width
                p['color']=color
                p['objLabelSize']=fontSize
                alreadyGot=True

        if alreadyGot == False:
            self.plotObjects.append({'x': [x], 'y': [y], 'RA': [RADeg], 'dec': [decDeg],
                                'tag': "compass", 'objLabels': [None], 'symbol': "compass",
                                'width': width, 'color': color,
                                'objLabelSize': fontSize, 'sizeArcSec': sizeArcSec})
        self.draw()


    def addScaleBar(self, location, sizeArcSec, color = "white", fontSize = 12, \
                        width = 20.0, label = None, style = "whiskers"):
        """Adds a scale bar to the ImagePlot at the given location ('N', 'NE', 'E', 'SE', 'S',
        'SW', 'W', or 'NW'). Note these aren't directions on the WCS coordinate grid, they are
        relative positions on the plot - so N is top centre, NE is top right, SW is bottom right etc..
        Alternatively, pixel coordinates (x, y) in the image can be given.

        Args:
            location (str or tuple): location in the plot where the scale bar is drawn:

                - string: N, NE, E, SE, S, SW, W or NW
                - tuple: (x, y)

            sizeArcSec (float): scale length to indicate on the plot in arc seconds
            color (str): any valid matplotlib color string
            fontSize (float): size of font used to label the scale bar, in points
            width (float): width of arrow used to mark scale
            label (str): overrides the displayed label if not None (if None, label is the angular size)
            style (str): either "whiskers" or "arrows"

        """

        # Work out where the scale bar is going in WCS coords from the relative location given
        if type(location) == str:
            cRADeg, cDecDeg=self.wcs.getCentreWCSCoords()
            RAMin, RAMax, decMin, decMax=self.wcs.getImageMinMaxWCSCoords()
            westPoint,eastPoint,southPoint,northPoint=astCoords.calcRADecSearchBox(cRADeg, cDecDeg, sizeArcSec/3600.0/2.0)
            sizeRADeg=eastPoint-westPoint
            sizeDecDeg=northPoint-southPoint
            xSizePix=(sizeArcSec/3600.0)/self.wcs.getXPixelSizeDeg()
            ySizePix=(sizeArcSec/3600.0)/self.wcs.getYPixelSizeDeg()
            X=self.data.shape[1]
            Y=self.data.shape[0]
            xBufferPix=0.6*ySizePix
            yBufferPix=0.05*Y
            cx, cy=self.wcs.wcs2pix(cRADeg, cDecDeg)
            foundLocation=False
            x=cy
            y=cx
            if self.wcs.isFlipped() == False:
                if location.find("N") != -1:
                    y=Y-1.5*yBufferPix
                    foundLocation=True
                if location.find("S") != -1:
                    y=yBufferPix
                    foundLocation=True
                if location.find("E") != -1:
                    x=xBufferPix
                    foundLocation=True
                if location.find("W") != -1:
                    x=X-xBufferPix
                    foundLocation=True
            else:
                if location.find("S") != -1:
                    y=Y-1.5*yBufferPix
                    foundLocation=True
                if location.find("N") != -1:
                    y=yBufferPix
                    foundLocation=True
                if location.find("W") != -1:
                    x=xBufferPix
                    foundLocation=True
                if location.find("E") != -1:
                    x=X-xBufferPix
                    foundLocation=True
            if foundLocation == False:
                raise Exception("didn't understand location string for scale bar (should be e.g. N, S, E, W).")
            RADeg, decDeg=self.wcs.pix2wcs(x, y)
        elif type(location) == tuple or type(location) == list:
            x, y=location
            RADeg, decDeg=self.wcs.pix2wcs(x, y)
        else:
            raise Exception("didn't understand location for scale bar - should be string or tuple.")

        alreadyGot=False
        for p in self.plotObjects:
            if p['tag'] == "scaleBar":
                p['x']=[x]
                p['y']=[y]
                p['RA']=[RADeg]
                p['dec']=[decDeg]
                p['tag']="scaleBar"
                p['objLabels']=[None]
                p['symbol']="scaleBar"
                p['sizeArcSec']=sizeArcSec
                p['width']=width
                p['color']=color
                p['objLabelSize']=fontSize
                p['scaleBarLabel']=label
                p['style']=style
                alreadyGot=True

        if alreadyGot == False:
            self.plotObjects.append({'x': [x], 'y': [y], 'RA': [RADeg], 'dec': [decDeg],
                                'tag': "scaleBar", 'objLabels': [None], 'symbol': "scaleBar",
                                'width': width, 'color': color,
                                'objLabelSize': fontSize, 'sizeArcSec': sizeArcSec,
                                'scaleBarLabel': label, 'style': style})
        self.draw()


    def calcWCSAxisLabels(self, axesLabels = "decimal"):
        """This function calculates the positions of coordinate labels for the RA and Dec axes of the
        ImagePlot. The tick steps are calculated automatically unless self.RATickSteps,
        self.decTickSteps are set to values other than "auto" (see __init__).

        The ImagePlot must be redrawn for changes to be applied.

        Args:
            axesLabels (str): either "sexagesimal" (for H:M:S, D:M:S), "decimal" (for decimal degrees),
                or None for no coordinate axes labels

        """

        # Label equinox on axes
        equinox=self.wcs.getEquinox()
        if equinox<1984:
            equinoxLabel="B"+str(int(equinox))
        else:
            equinoxLabel="J"+str(int(equinox))

        self.axesLabels=axesLabels

        ticsDict=self.getTickSteps()

        # Manual override - note: no minor tick marks anymore, but may want to bring them back
        if self.RATickSteps != "auto":
            ticsDict['major']['RA']=self.RATickSteps
        if self.decTickSteps != "auto":
            ticsDict['major']['dec']=self.decTickSteps

        RALocs=[]
        decLocs=[]
        RALabels=[]
        decLabels=[]
        key="major"
        #for key in ticsDict.keys(): # key is major or minor
        if self.axesLabels == "sexagesimal":
            self.RAAxisLabel="R.A. ("+equinoxLabel+")"
            self.decAxisLabel="Dec. ("+equinoxLabel+")"
            RADegStep=ticsDict[key]['RA']['deg']
            decDegStep=ticsDict[key]['dec']['deg']
        elif self.axesLabels == "decimal":
            self.RAAxisLabel="R.A. Degrees ("+equinoxLabel+")"
            self.decAxisLabel="Dec. Degrees ("+equinoxLabel+")"
            RADegStep=ticsDict[key]['RA']
            decDegStep=ticsDict[key]['dec']
        else:
            raise Exception("axesLabels must be either 'sexagesimal' or 'decimal'")

        # xArray=np.arange(0, self.data.shape[1], 1)
        xArray=np.linspace(0, self.data.shape[1], self.data.shape[1]*2)
        yArray=np.arange(0, self.data.shape[0], 1)
        xWCS=self.wcs.pix2wcs(xArray, np.zeros(xArray.shape[0], dtype=float))
        yWCS=self.wcs.pix2wcs(np.zeros(yArray.shape[0], dtype=float), yArray)
        xWCS=np.array(xWCS)
        yWCS=np.array(yWCS)
        ras=xWCS[:,0]
        # ras[ras < 0]=ras[ras < 0]+360
        decs=yWCS[:,1]
        RAEdges=np.array([ras[0], ras[-1]])
        # RAMin=RAEdges.min()
        # RAMax=RAEdges.max()
        RAMin=ras.min()
        RAMax=ras.max()
        decMin=decs.min()
        decMax=decs.max()

        # Work out if wrapped around
        midRAPix, midDecPix=self.wcs.wcs2pix((RAEdges[1]+RAEdges[0])/2.0, (decMax+decMin)/2.0)
        if midRAPix < 0 or midRAPix > self.wcs.header['NAXIS1']:
            wrappedRA=True
        else:
            wrappedRA=False

        # Note RA, dec work in opposite sense below because E at left
        if ras[1] < ras[0]:
            self.flipXAxis=False
            ra2x=interpolate.interp1d(ras[::-1], xArray[::-1], kind='linear')
        else:
            self.flipXAxis=True
            ra2x=interpolate.interp1d(ras, xArray, kind='linear')
        if decs[1] < decs[0]:
            self.flipYAxis=True
            dec2y=interpolate.interp1d(decs[::-1], yArray[::-1], kind='linear')
        else:
            self.flipYAxis=False
            dec2y=interpolate.interp1d(decs, yArray, kind='linear')

        if wrappedRA == False:
            RAPlotMin=RADegStep*math.modf(RAMin/RADegStep)[1]
            RAPlotMax=RADegStep*math.modf(RAMax/RADegStep)[1]
            if RAPlotMin < RAMin:
                RAPlotMin=RAPlotMin+RADegStep
            if RAPlotMax >= RAMax:
                RAPlotMax=RAPlotMax-RADegStep
            RADegs=np.arange(RAPlotMin, RAPlotMax+0.0001, RADegStep)
        else:
            RAPlotMin=RADegStep*math.modf(RAMin/RADegStep)[1]
            RAPlotMax=RADegStep*math.modf(RAMax/RADegStep)[1]
            if RAPlotMin > RAMin:
                RAPlotMin=RAPlotMin-RADegStep
            if RAPlotMax <= RAMax:
                RAPlotMax=RAPlotMax+RADegStep
            for i in range(ras.shape[0]):
                if ras[i] >= RAMax and ras[i] <= 360.0:
                    ras[i]=ras[i]-360.0
            if ras[1] < ras[0]:
                ra2x=interpolate.interp1d(ras[::-1], xArray[::-1], kind='linear')
            else:
                ra2x=interpolate.interp1d(ras, xArray, kind='linear')
            RADegs=np.arange(RAPlotMin, RAPlotMax-360.0-0.0001, -RADegStep)

        # Full RA range?
        if round(self.wcs.getXPixelSizeDeg()*self.wcs.header['NAXIS1']) == 360:
            RADegs=np.arange(0, 360, RADegStep)

        decPlotMin=decDegStep*math.modf(decMin/decDegStep)[1]
        decPlotMax=decDegStep*math.modf(decMax/decDegStep)[1]
        if decPlotMin < decMin:
            decPlotMin=decPlotMin+decDegStep
        if decPlotMax >= decMax:
            decPlotMax=decPlotMax-decDegStep
        decDegs=np.arange(decPlotMin, decPlotMax+0.0001, decDegStep)

        if key == "major":
            if axesLabels == "sexagesimal":
                for r in RADegs:
                    if r < 0:
                        r=r+360.0
                    h, m, s=astCoords.decimal2hms(r, ":").split(":")
                    hInt=int(round(float(h)))
                    if ticsDict[key]['RA']['unit'] == 'h' and (60.0-float(m)) < 0.01: # Check for rounding error
                        hInt=hInt+1
                    if hInt < 10:
                        hString="0"+str(hInt)
                    else:
                        hString=str(hInt)
                    mInt=int(round(float(m)))
                    if ticsDict[key]['RA']['unit'] == 'm' and (60.0-float(s)) < 0.01: # Check for rounding error
                        mInt=mInt+1
                    if mInt < 10:
                        mString="0"+str(mInt)
                    else:
                        mString=str(mInt)
                    sInt=int(round(float(s)))
                    if sInt < 10:
                        sString="0"+str(sInt)
                    else:
                        sString=str(sInt)
                    if ticsDict[key]['RA']['unit'] == 'h':
                        rString=hString+"$^{\sf{h}}$"
                    elif ticsDict[key]['RA']['unit'] == 'm':
                        rString=hString+"$^{\sf{h}}$"+mString+"$^{\sf{m}}$"
                    else:
                        rString=hString+"$^{\sf{h}}$"+mString+"$^{\sf{m}}$"+sString+"$^{\sf{s}}$"
                    RALabels.append(rString)
                for D in decDegs:
                    d, m, s=astCoords.decimal2dms(D, ":").split(":")
                    dInt=int(round(float(d)))
                    if ticsDict[key]['dec']['unit'] == 'd' and (60.0-float(m)) < 0.01: # Check for rounding error
                        dInt=dInt+1
                    if dInt < 10 and dInt >= 0 and D > 0:
                        dString="+0"+str(dInt)
                    elif dInt > -10 and dInt <= 0 and D < 0:
                        dString="-0"+str(abs(dInt))
                    elif dInt >= 10:
                        dString="+"+str(dInt)
                    else:
                        dString=str(dInt)
                    mInt=int(round(float(m)))
                    if ticsDict[key]['dec']['unit'] == 'm' and (60.0-float(s)) < 0.01: # Check for rounding error
                        mInt=mInt+1
                    if mInt < 10:
                        mString="0"+str(mInt)
                    else:
                        mString=str(mInt)
                    sInt=int(round(float(s)))
                    if sInt < 10:
                        sString="0"+str(sInt)
                    else:
                        sString=str(sInt)
                    if ticsDict[key]['dec']['unit'] == 'd':
                        dString=dString+DEG
                    elif ticsDict[key]['dec']['unit'] == 'm':
                        dString=dString+DEG+mString+PRIME
                    else:
                        dString=dString+DEG+mString+PRIME+sString+DOUBLE_PRIME
                    decLabels.append(dString)
            elif axesLabels == "decimal":

                if wrappedRA == False:
                    RALabels=RALabels+RADegs.tolist()
                else:
                    nonNegativeLabels=[]
                    for r in RADegs:
                        if r < 0:
                            r=r+360.0
                        nonNegativeLabels.append(r)
                    RALabels=RALabels+nonNegativeLabels
                decLabels=decLabels+decDegs.tolist()

                # Format RALabels, decLabels to same number of d.p.
                dpNumRA=len(str(ticsDict['major']['RA']).split(".")[-1])
                dpNumDec=len(str(ticsDict['major']['dec']).split(".")[-1])
                for i in range(len(RALabels)):
                    fString="%."+str(dpNumRA)+"f"
                    RALabels[i]=fString % (RALabels[i])
                for i in range(len(decLabels)):
                    fString="%."+str(dpNumDec)+"f"
                    decLabels[i]=fString % (decLabels[i])

        if key == 'minor':
            RALabels=RALabels+RADegs.shape[0]*['']
            decLabels=decLabels+decDegs.shape[0]*['']

        decLocs=decLocs+dec2y(decDegs).tolist()

        # This is for avoiding 0h <-> 12h degeneracy that sometimes happens
        RALocs=[]
        RALocs=RALocs+ra2x(RADegs).tolist()
        uniqRALocs=[]
        uniqRALabels=[]
        for i in range(len(RALocs)):
            if RALocs[i] not in uniqRALocs:
                uniqRALocs.append(RALocs[i])
                uniqRALabels.append(RALabels[i])
        RALocs=uniqRALocs
        RALabels=uniqRALabels

        self.ticsRA=[RALocs, RALabels]
        self.ticsDec=[decLocs, decLabels]


    def save(self, fileName):
        """Saves the ImagePlot in any format that matplotlib can understand, as determined from the
        fileName extension.

        Args:
            fileName (str): path where plot will be written

        """

        pylab.draw()
        pylab.savefig(fileName)


    def getTickSteps(self):
        """Chooses the appropriate WCS coordinate tick steps for the plot based on its size.
        Whether the ticks are decimal or sexagesimal is set by self.axesLabels.

        Note:
            Minor ticks not used at the moment.

        Returns:
            dict: tick step sizes for major, minor plot ticks, in format ``{'major', 'minor'}``

        """

        # Aim for 5 major tick marks on a plot
        xArray=np.arange(0, self.data.shape[1], 1)
        yArray=np.arange(0, self.data.shape[0], 1)
        xWCS=self.wcs.pix2wcs(xArray, np.zeros(xArray.shape[0], dtype=float))
        yWCS=self.wcs.pix2wcs(np.zeros(yArray.shape[0], dtype=float), yArray)
        xWCS=np.array(xWCS)
        yWCS=np.array(yWCS)
        ras=xWCS[:,0]
        decs=yWCS[:,1]
        RAEdges=np.array([ras[0], ras[-1]])
        RAMin=RAEdges.min()
        RAMax=RAEdges.max()
        decMin=decs.min()
        decMax=decs.max()

        # Work out if wrapped around
        midRAPix, midDecPix=self.wcs.wcs2pix((RAEdges[1]+RAEdges[0])/2.0, (decMax+decMin)/2.0)
        if midRAPix < 0 or midRAPix > self.wcs.header['NAXIS1']:
            wrappedRA=True
        else:
            wrappedRA=False
        if wrappedRA == False:
            RAWidthDeg=RAMax-RAMin
        else:
            RAWidthDeg=(360.0-RAMax)+RAMin
        decHeightDeg=decMax-decMin

        ticsDict={}
        ticsDict['major']={}
        ticsDict['minor']={}
        if self.axesLabels == "sexagesimal":

            matchIndex = 0
            for i in range(len(RA_TICK_STEPS)):
                if RAWidthDeg/2.5 > RA_TICK_STEPS[i]['deg']:
                    matchIndex = i

            ticsDict['major']['RA']=RA_TICK_STEPS[matchIndex]
            ticsDict['minor']['RA']=RA_TICK_STEPS[matchIndex-1]

            matchIndex = 0
            for i in range(len(DEC_TICK_STEPS)):
                if decHeightDeg/2.5 > DEC_TICK_STEPS[i]['deg']:
                    matchIndex = i

            ticsDict['major']['dec']=DEC_TICK_STEPS[matchIndex]
            ticsDict['minor']['dec']=DEC_TICK_STEPS[matchIndex-1]

            return ticsDict

        elif self.axesLabels == "decimal":

            matchIndex = 0
            for i in range(len(DECIMAL_TICK_STEPS)):
                if RAWidthDeg/2.5 > DECIMAL_TICK_STEPS[i]:
                    matchIndex = i

            ticsDict['major']['RA']=DECIMAL_TICK_STEPS[matchIndex]
            ticsDict['minor']['RA']=DECIMAL_TICK_STEPS[matchIndex-1]

            matchIndex = 0
            for i in range(len(DECIMAL_TICK_STEPS)):
                if decHeightDeg/2.5 > DECIMAL_TICK_STEPS[i]:
                    matchIndex = i

            ticsDict['major']['dec']=DECIMAL_TICK_STEPS[matchIndex]
            ticsDict['minor']['dec']=DECIMAL_TICK_STEPS[matchIndex-1]

            return ticsDict

        else:
            raise Exception("axesLabels must be either 'sexagesimal' or 'decimal'")
