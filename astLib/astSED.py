"""module for performing calculations on Spectral Energy Distributions (SEDs)

(c) 2007-2013 Matt Hilton 

U{http://astlib.sourceforge.net}

This module provides classes for manipulating SEDs, in particular the Bruzual & Charlot 2003, Maraston
2005, and Percival et al 2009 stellar population synthesis models are currently supported. Functions are 
provided for calculating the evolution of colours and magnitudes in these models with redshift etc., and 
for fitting broadband photometry using these models.

@var VEGA: The SED of Vega, used for calculation of magnitudes on the Vega system.
@type VEGA: L{SED} object
@var AB: Flat spectrum SED, used for calculation of magnitudes on the AB system.
@type AB: L{SED} object
@var SOL: The SED of the Sun.
@type SOL: L{SED} object

"""

#------------------------------------------------------------------------------------------------------------
import sys
import numpy
import math
import operator
try:
    from scipy import interpolate
    from scipy import ndimage
    from scipy import optimize
except:
    print("WARNING: astSED: failed to import scipy modules - some functions will not work.")
import astLib
from astLib import astCalc
import os
try:
    import matplotlib
    from matplotlib import pylab
    matplotlib.interactive(False)
except:
    print("WARNING: astSED: failed to import matplotlib - some functions will not work.")
import glob

#------------------------------------------------------------------------------------------------------------
class Passband:
    """This class describes a filter transmission curve. Passband objects are created by loading data from
    from text files containing wavelength in angstroms in the first column, relative transmission efficiency
    in the second column (whitespace delimited). For example, to create a Passband object for the 2MASS J 
    filter:
    
    passband=astSED.Passband("J_2MASS.res")
    
    where "J_2MASS.res" is a file in the current working directory that describes the filter.
    
    Wavelength units can be specified as 'angstroms', 'nanometres' or 'microns'; if either of the latter,
    they will be converted to angstroms.
    
    """
    def __init__(self, fileName, normalise = True, inputUnits = 'angstroms', wavelengthColumn = 0, transmissionColumn = 1):
        
        inFile=open(fileName, "r")
        lines=inFile.readlines()
        
        wavelength=[]
        transmission=[]
        for line in lines:
            
            if line[0] != "#" and len(line) > 3:
    
                bits=line.split()
                transmission.append(float(bits[transmissionColumn]))
                wavelength.append(float(bits[wavelengthColumn]))
            
        self.wavelength=numpy.array(wavelength)
        self.transmission=numpy.array(transmission)
        
        if inputUnits == 'angstroms':
            pass
        elif inputUnits == 'nanometres':
            self.wavelength=self.wavelength*10.0
        elif inputUnits == 'microns':
            self.wavelength=self.wavelength*10000.0
        elif inputUnits == 'mm':
            self.wavelength=self.wavelength*1e7
        elif inputUnits == 'GHz':
            self.wavelength=3e8/(self.wavelength*1e9)
            self.wavelength=self.wavelength*1e10
        else:
            raise Exception("didn't understand passband input units")
    
        # Sort into ascending order of wavelength otherwise normalisation will be wrong
        merged=numpy.array([self.wavelength, self.transmission]).transpose()
        sortedMerged=numpy.array(sorted(merged, key=operator.itemgetter(0)))
        self.wavelength=sortedMerged[:, 0]
        self.transmission=sortedMerged[:, 1]
        
        if normalise == True:
            self.transmission=self.transmission/numpy.trapz(self.transmission, self.wavelength)
        
        # Store a ready-to-go interpolation object to speed calculation of fluxes up
        self.interpolator=interpolate.interp1d(self.wavelength, self.transmission, kind='linear')

    def asList(self):
        """Returns a two dimensional list of [wavelength, transmission], suitable for plotting by gnuplot.
        
        @rtype: list
        @return: list in format [wavelength, transmission]
        
        """
        
        listData=[]
        for l, f in zip(self.wavelength, self.transmission):
            listData.append([l, f])
        
        return listData
    
    def rescale(self, maxTransmission):
        """Rescales the passband so that maximum value of the transmission is equal to maxTransmission.
        Useful for plotting.
        
        @type maxTransmission: float
        @param maxTransmission: maximum value of rescaled transmission curve
        
        """
        
        self.transmission=self.transmission*(maxTransmission/self.transmission.max())

    def plot(self, xmin = 'min', xmax = 'max', maxTransmission = None):
        """Plots the passband, rescaling the maximum of the tranmission curve to maxTransmission if
        required.
        
        @type xmin: float or 'min'
        @param xmin: minimum of the wavelength range of the plot
        @type xmax: float or 'max'
        @param xmax: maximum of the wavelength range of the plot
        @type maxTransmission: float
        @param maxTransmission: maximum value of rescaled transmission curve
        
        """
        
        if maxTransmission != None:
            self.rescale(maxTransmission)
        
        pylab.matplotlib.interactive(True)
        pylab.plot(self.wavelength, self.transmission)
        
        if xmin == 'min':
            xmin=self.wavelength.min()
        if xmax == 'max':
            xmax=self.wavelength.max()
            
        pylab.xlim(xmin, xmax)
        pylab.xlabel("Wavelength")
        pylab.ylabel("Relative Flux")

    def effectiveWavelength(self):
        """Calculates effective wavelength for the passband. This is the same as equation (3) of
        Carter et al. 2009.
        
        @rtype: float
        @return: effective wavelength of the passband, in Angstroms
        
        """
        
        a=numpy.trapz(self.transmission*self.wavelength, self.wavelength)
        b=numpy.trapz(self.transmission/self.wavelength, self.wavelength)
        effWavelength=numpy.sqrt(a/b)
        
        return effWavelength

#------------------------------------------------------------------------------------------------------------
class TopHatPassband(Passband):
    """This class generates a passband with a top hat response between the given wavelengths.
    
    """
    
    def __init__(self, wavelengthMin, wavelengthMax, normalise = True):
        """Generates a passband object with top hat response between wavelengthMin, wavelengthMax.
        Units are assumed to be Angstroms.
        
        @type wavelengthMin: float
        @param wavelengthMin: minimum of the wavelength range of the passband
        @type wavelengthMax: float
        @param wavelengthMax: maximum of the wavelength range of the passband
        @type normalise: bool
        @param normalise: if True, scale such that total area under the passband over the wavelength 
        range is 1.
        
        """
        
        self.wavelength=numpy.arange(wavelengthMin, wavelengthMax+10, 10, dtype = float)
        self.transmission=numpy.ones(self.wavelength.shape, dtype = float)
        
        if normalise == True:
            self.transmission=self.transmission/numpy.trapz(self.transmission, self.wavelength)
        
        # Store a ready-to-go interpolation object to speed calculation of fluxes up
        self.interpolator=interpolate.interp1d(self.wavelength, self.transmission, kind='linear')
        
    
#------------------------------------------------------------------------------------------------------------
class SED:
    """This class describes a Spectral Energy Distribution (SED).
     
    To create a SED object, lists (or numpy arrays) of wavelength and relative flux must be provided. The SED
    can optionally be redshifted. The wavelength units of SEDs are assumed to be Angstroms - flux 
    calculations using Passband and SED objects specified with different wavelength units will be incorrect.
    
    The L{StellarPopulation} class (and derivatives) can be used to extract SEDs for specified ages from e.g.
    the Bruzual & Charlot 2003 or Maraston 2005 models.
    
    """
    
    def __init__(self, wavelength = [], flux = [], z = 0.0, ageGyr = None, normalise = False, label = None):
        
        # We keep a copy of the wavelength, flux at z = 0, as it's more robust to copy that
        # to self.wavelength, flux and redshift it, rather than repeatedly redshifting the same
        # arrays back and forth
        self.z0wavelength=numpy.array(wavelength)
        self.z0flux=numpy.array(flux)
        self.wavelength=numpy.array(wavelength)
        self.flux=numpy.array(flux)
        self.z=z
        self.label=label    # plain text label, handy for using in photo-z codes
        
        # Store the intrinsic (i.e. unextincted) flux in case we change extinction
        self.EBMinusV=0.0
        self.intrinsic_z0flux=numpy.array(flux)
        
        if normalise == True:
            self.normalise()
            
        if z != 0.0:
            self.redshift(z)

        self.ageGyr=ageGyr


    def copy(self):
        """Copies the SED, returning a new SED object
        
        @rtype: L{SED} object
        @return: SED
        
        """
        
        newSED=SED(wavelength = self.z0wavelength, flux = self.z0flux, z = self.z, ageGyr = self.ageGyr, 
                   normalise = False, label = self.label)
        
        return newSED
        
        
    def loadFromFile(self, fileName):
        """Loads SED from a white space delimited file in the format wavelength, flux. Lines beginning with
        # are ignored.
        
        @type fileName: string
        @param fileName: path to file containing wavelength, flux data
        
        """
        
        inFile=open(fileName, "r")
        lines=inFile.readlines()
        inFile.close()
        wavelength=[]
        flux=[]
        wholeLines=[]
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                bits=line.split()
                wavelength.append(float(bits[0]))
                flux.append(float(bits[1]))
        
        # Sort SED so wavelength is in ascending order
        if wavelength[0] > wavelength[-1]:
            wavelength.reverse()
            flux.reverse()
        
        self.z0wavelength=numpy.array(wavelength)
        self.z0flux=numpy.array(flux)
        self.wavelength=numpy.array(wavelength)
        self.flux=numpy.array(flux)

    def writeToFile(self, fileName):
        """Writes SED to a white space delimited file in the format wavelength, flux.
        
        @type fileName: string
        @param fileName: path to file
        
        """
        
        outFile=open(fileName, "w")
        for l, f in zip(self.wavelength, self.flux):
            outFile.write(str(l)+" "+str(f)+"\n")
        outFile.close()
    
    def asList(self):
        """Returns a two dimensional list of [wavelength, flux], suitable for plotting by gnuplot.
        
        @rtype: list
        @return: list in format [wavelength, flux]
        
        """
        
        listData=[]
        for l, f in zip(self.wavelength, self.flux):
            listData.append([l, f])
        
        return listData
        
    def plot(self, xmin = 'min', xmax = 'max'):
        """Produces a simple (wavelength, flux) plot of the SED.
        
        @type xmin: float or 'min'
        @param xmin: minimum of the wavelength range of the plot
        @type xmax: float or 'max'
        @param xmax: maximum of the wavelength range of the plot
        
        """
        
        pylab.matplotlib.interactive(True)
        pylab.plot(self.wavelength, self.flux)
        
        if xmin == 'min':
            xmin=self.wavelength.min()
        if xmax == 'max':
            xmax=self.wavelength.max()
        
        # Sensible y scale
        plotMask=numpy.logical_and(numpy.greater(self.wavelength, xmin), numpy.less(self.wavelength, xmax))
        plotMax=self.flux[plotMask].max()
        pylab.ylim(0, plotMax*1.1)
        pylab.xlim(xmin, xmax)
        pylab.xlabel("Wavelength")
        pylab.ylabel("Relative Flux")
    
    def integrate(self, wavelengthMin = 'min', wavelengthMax = 'max'):
        """Calculates flux in SED within given wavelength range.
        
        @type wavelengthMin: float or 'min'
        @param wavelengthMin: minimum of the wavelength range
        @type wavelengthMax: float or 'max'
        @param wavelengthMax: maximum of the wavelength range
        @rtype: float
        @return: relative flux
        
        """

        if wavelengthMin == 'min':
            wavelengthMin=self.wavelength.min()
        if wavelengthMax == 'max':
            wavelengthMax=self.wavelength.max()
        
        mask=numpy.logical_and(numpy.greater(self.wavelength, wavelengthMin), \
                               numpy.less(self.wavelength, wavelengthMax))
        flux=numpy.trapz(self.flux[mask], self.wavelength[mask])
        
        return flux
        
    def smooth(self, smoothPix):
        """Smooths SED.flux with a uniform (boxcar) filter of width smoothPix. Cannot be undone.
        
        @type smoothPix: int
        @param smoothPix: size of uniform filter applied to SED, in pixels
        
        """
        smoothed=ndimage.uniform_filter1d(self.flux, smoothPix)
        self.flux=smoothed
    
    def redshift(self, z):
        """Redshifts the SED to redshift z.
        
        @type z: float
        @param z: redshift
        
        """
        
        # We have to conserve energy so the area under the redshifted SED has to be equal to
        # the area under the unredshifted SED, otherwise magnitude calculations will be wrong
        # when comparing SEDs at different zs
        self.wavelength=numpy.zeros(self.z0wavelength.shape[0])
        self.flux=numpy.zeros(self.z0flux.shape[0])
        self.wavelength=self.wavelength+self.z0wavelength
        self.flux=self.flux+self.z0flux
        
        z0TotalFlux=numpy.trapz(self.z0wavelength, self.z0flux)
        self.wavelength=self.wavelength*(1.0+z)
        zTotalFlux=numpy.trapz(self.wavelength, self.flux)
        self.flux=self.flux*(z0TotalFlux/zTotalFlux)
        self.z=z
        
    def normalise(self, minWavelength = 'min', maxWavelength = 'max'):
        """Normalises the SED such that the area under the specified wavelength range is equal to 1.
        
        @type minWavelength: float or 'min'
        @param minWavelength: minimum wavelength of range over which to normalise SED
        @type maxWavelength: float or 'max'
        @param maxWavelength: maximum wavelength of range over which to normalise SED
        
        """
        if minWavelength == 'min':
            minWavelength=self.wavelength.min()
        if maxWavelength == 'max':
            maxWavelength=self.wavelength.max()
            
        lowCut=numpy.greater(self.wavelength, minWavelength)
        highCut=numpy.less(self.wavelength, maxWavelength)
        totalCut=numpy.logical_and(lowCut, highCut)
        sedFluxSlice=self.flux[totalCut]
        sedWavelengthSlice=self.wavelength[totalCut]
        
        self.flux=self.flux/numpy.trapz(abs(sedFluxSlice), sedWavelengthSlice)#self.wavelength)

    def normaliseToMag(self, ABMag, passband):
        """Normalises the SED to match the flux equivalent to the given AB magnitude in the given passband.
        
        @type ABMag: float
        @param ABMag: AB magnitude to which the SED is to be normalised at the given passband
        @type passband: an L{Passband} object
        @param passband: passband at which normalisation to AB magnitude is calculated
        
        """
        
        magFlux=mag2Flux(ABMag, 0.0, passband)
        sedFlux=self.calcFlux(passband)
        norm=magFlux[0]/sedFlux
        self.flux=self.flux*norm
        self.z0flux=self.z0flux*norm
        
    def matchFlux(self, matchSED, minWavelength, maxWavelength):
        """Matches the flux in the wavelength range given by minWavelength, maxWavelength to the
        flux in the same region in matchSED. Useful for plotting purposes.
        
        @type matchSED: L{SED} object
        @param matchSED: SED to match flux to
        @type minWavelength: float
        @param minWavelength: minimum of range in which to match flux of current SED to matchSED
        @type maxWavelength: float
        @param maxWavelength: maximum of range in which to match flux of current SED to matchSED
        
        """
        
        interpMatch=interpolate.interp1d(matchSED.wavelength, matchSED.flux, kind='linear')
        interpSelf=interpolate.interp1d(self.wavelength, self.flux, kind='linear')
        
        wavelengthRange=numpy.arange(minWavelength, maxWavelength, 5.0)
        
        matchFlux=numpy.trapz(interpMatch(wavelengthRange), wavelengthRange)
        selfFlux=numpy.trapz(interpSelf(wavelengthRange), wavelengthRange)
        
        self.flux=self.flux*(matchFlux/selfFlux)

        
    def calcFlux(self, passband):
        """Calculates flux in the given passband.
        
        @type passband: L{Passband} object
        @param passband: filter passband through which to calculate the flux from the SED
        @rtype: float
        @return: flux
        
        """
        lowCut=numpy.greater(self.wavelength, passband.wavelength.min())
        highCut=numpy.less(self.wavelength, passband.wavelength.max())
        totalCut=numpy.logical_and(lowCut, highCut)
        sedFluxSlice=self.flux[totalCut]
        sedWavelengthSlice=self.wavelength[totalCut]
    
        # Use linear interpolation to rebin the passband to the same dimensions as the 
        # part of the SED we're interested in
        sedInBand=passband.interpolator(sedWavelengthSlice)*sedFluxSlice   
        totalFlux=numpy.trapz(sedInBand*sedWavelengthSlice, sedWavelengthSlice)        
        totalFlux=totalFlux/numpy.trapz(passband.interpolator(sedWavelengthSlice)\
                            *sedWavelengthSlice, sedWavelengthSlice)
                            
        return totalFlux      
    
    def calcMag(self, passband, addDistanceModulus = True, magType = "Vega"):
        """Calculates magnitude in the given passband. If addDistanceModulus == True,
        then the distance modulus (5.0*log10*(dl*1e5), where dl is the luminosity distance
        in Mpc at the redshift of the L{SED}) is added.
               
        @type passband: L{Passband} object
        @param passband: filter passband through which to calculate the magnitude from the SED
        @type addDistanceModulus: bool
        @param addDistanceModulus: if True, adds 5.0*log10*(dl*1e5) to the mag returned, where
                                   dl is the luminosity distance (Mpc) corresponding to the SED z
        @type magType: string
        @param magType: either "Vega" or "AB"
        @rtype: float
        @return: magnitude through the given passband on the specified magnitude system
        
        """
        f1=self.calcFlux(passband)
        if magType == "Vega":
            f2=VEGA.calcFlux(passband)
        elif magType == "AB":
            f2=AB.calcFlux(passband)
        
        mag=-2.5*math.log10(f1/f2)
        if magType == "Vega":
            mag=mag+0.026               # Add 0.026 because Vega has V=0.026 (e.g. Bohlin & Gilliland 2004)
                    
        if self.z > 0.0 and addDistanceModulus == True:
            appMag=5.0*math.log10(astCalc.dl(self.z)*1e5)+mag
        else:
            appMag=mag
        
        return appMag
    
    def calcColour(self, passband1, passband2, magType = "Vega"):
        """Calculates the colour passband1-passband2.
        
        @type passband1: L{Passband} object
        @param passband1: filter passband through which to calculate the first magnitude
        @type passband2: L{Passband} object
        @param passband1: filter passband through which to calculate the second magnitude
        @type magType: string
        @param magType: either "Vega" or "AB"
        @rtype: float
        @return: colour defined by passband1 - passband2 on the specified magnitude system

        """
        mag1=self.calcMag(passband1, magType = magType, addDistanceModulus = True)
        mag2=self.calcMag(passband2, magType = magType, addDistanceModulus = True)
        
        colour=mag1-mag2
        return colour
    
    def getSEDDict(self, passbands):
        """This is a convenience function for pulling out fluxes from a SED for a given set of passbands
        in the same format as made by L{mags2SEDDict} - designed to make fitting code simpler.
        
        @type passbands: list of L{Passband} objects
        @param passbands: list of passbands through which fluxes will be calculated
        
        """
        
        flux=[]
        wavelength=[]
        for p in passbands:
            flux.append(self.calcFlux(p))
            wavelength.append(p.effectiveWavelength())
            
        SEDDict={}
        SEDDict['flux']=numpy.array(flux)
        SEDDict['wavelength']=numpy.array(wavelength)
        
        return SEDDict
    
    def extinctionCalzetti(self, EBMinusV):
        """Applies the Calzetti et al. 2000 (ApJ, 533, 682) extinction law to the SED with the given
        E(B-V) amount of extinction. R_v' = 4.05 is assumed (see equation (5) of Calzetti et al.).
        
        @type EBMinusV: float
        @param EBMinusV: extinction E(B-V), in magnitudes
        
        """
        
        self.EBMinusV=EBMinusV
        
        # All done in rest frame
        self.z0flux=self.intrinsic_z0flux
        
        # Allow us to set EBMinusV == 0 to turn extinction off
        if EBMinusV > 0:
            # Note that EBMinusV is assumed to be Es as in equations (2) - (5)
            # Note here wavelength units have to be microns for constants to make sense
            RvPrime=4.05    # equation (5) of Calzetti et al. 2000
            shortWavelengthMask=numpy.logical_and(numpy.greater_equal(self.z0wavelength, 1200), \
                                                 numpy.less(self.z0wavelength, 6300))
            longWavelengthMask=numpy.logical_and(numpy.greater_equal(self.z0wavelength, 6300), \
                                                numpy.less_equal(self.z0wavelength, 22000))
            wavelengthMicrons=numpy.array(self.z0wavelength/10000.0, dtype=numpy.float64)
            kPrime=numpy.zeros(self.z0wavelength.shape[0], dtype=numpy.float64)
            kPrimeLong=(2.659*(-1.857 \
                                +1.040/wavelengthMicrons \
                               ))+RvPrime
            kPrimeShort=(2.659*(-2.156 \
                                +1.509/wavelengthMicrons \
                                -0.198/wavelengthMicrons**2 \
                                +0.011/wavelengthMicrons**3 \
                               ))+RvPrime
            kPrime[longWavelengthMask]=kPrimeLong[longWavelengthMask]
            kPrime[shortWavelengthMask]=kPrimeShort[shortWavelengthMask]

            # Here we extrapolate kPrime in similar way to what HYPERZ does
            # Short wavelengths
            try:
                interpolator=interpolate.interp1d(self.z0wavelength, kPrimeShort, kind='linear')
                slope=(interpolator(1100.0)-interpolator(1200.0))/(1100.0-1200.0)
                intercept=interpolator(1200.0)-(slope*1200.0)
                mask=numpy.less(self.z0wavelength, 1200.0)
                kPrime[mask]=slope*self.z0wavelength[mask]+intercept
                
                # Long wavelengths
                interpolator=interpolate.interp1d(self.z0wavelength, kPrimeLong, kind='linear')
                slope=(interpolator(21900.0)-interpolator(22000.0))/(21900.0-22000.0)
                intercept=interpolator(21900.0)-(slope*21900.0)
                mask=numpy.greater(self.z0wavelength, 22000.0)
                kPrime[mask]=slope*self.z0wavelength[mask]+intercept
            except:
                raise Exception("This SED has a wavelength range that doesn't cover ~1200-22000 Angstroms")
                            
            # Never let go negative
            kPrime[numpy.less_equal(kPrime, 0.0)]=1e-6
                
            reddening=numpy.power(10, 0.4*EBMinusV*kPrime)
            self.z0flux=self.z0flux/reddening

        self.redshift(self.z)
        
#------------------------------------------------------------------------------------------------------------
class VegaSED(SED):
    """This class stores the SED of Vega, used for calculation of magnitudes on the Vega system.
    
    The Vega SED used is taken from Bohlin 2007 (http://adsabs.harvard.edu/abs/2007ASPC..364..315B), and is
    available from the STScI CALSPEC library (http://www.stsci.edu/hst/observatory/cdbs/calspec.html).
    
    """
    
    def __init__(self, normalise = False):
        
        VEGA_SED_PATH=astLib.__path__[0]+os.path.sep+"data"+os.path.sep+"bohlin2006_Vega.sed" # from HST CALSPEC

        inFile=open(VEGA_SED_PATH, "r")
        lines=inFile.readlines()
        
        wavelength=[]
        flux=[]
        for line in lines:
            
            if line[0] != "#" and len(line) > 3:
        
                bits=line.split()
                flux.append(float(bits[1]))
                wavelength.append(float(bits[0]))
        
        self.wavelength=numpy.array(wavelength)
        self.flux=numpy.array(flux, dtype=numpy.float64)
        
        # We may want to redshift reference SEDs to calculate rest-frame colors from SEDs at different zs
        self.z0wavelength=numpy.array(wavelength)
        self.z0flux=numpy.array(flux, dtype=numpy.float64)
        self.z=0.0
        
        #if normalise == True:
            #self.flux=self.flux/numpy.trapz(self.flux, self.wavelength)
            #self.z0flux=self.z0flux/numpy.trapz(self.z0flux, self.z0wavelength)
        
#------------------------------------------------------------------------------------------------------------
class StellarPopulation:
    """This class describes a stellar population model, either a Simple Stellar Population (SSP) or a
    Composite Stellar Population (CSP), such as the models of Bruzual & Charlot 2003 or Maraston 2005.
    
    The constructor for this class can be used for generic SSPs or CSPs stored in white space delimited text
    files, containing columns for age, wavelength, and flux. Columns are counted from 0 ... n. Lines starting
    with # are ignored.
    
    The classes L{M05Model} (for Maraston 2005 models), L{BC03Model} (for Bruzual & Charlot 2003 models), and
    L{P09Model} (for Percival et al. 2009 models) are derived from this class. The only difference between 
    them is the code used to load in the model data.
   
    """
    def __init__(self, fileName, ageColumn = 0, wavelengthColumn = 1, fluxColumn = 2):
       
        inFile=open(fileName, "r")
        lines=inFile.readlines()
        inFile.close()

        self.fileName=fileName

        # Extract a list of model ages and valid wavelengths from the file
        self.ages=[]
        self.wavelengths=[]
        for line in lines:
            if line[0] !="#" and len(line) > 3:
                bits=line.split()
                age=float(bits[ageColumn])
                wavelength=float(bits[wavelengthColumn])
                if age not in self.ages:
                    self.ages.append(age)
                if wavelength not in self.wavelengths:
                    self.wavelengths.append(wavelength)
        
        # Construct a grid of flux - rows correspond to each wavelength, columns to age
        self.fluxGrid=numpy.zeros([len(self.ages), len(self.wavelengths)])
        for line in lines:
            if line[0] !="#" and len(line) > 3:
                bits=line.split()
                sedAge=float(bits[ageColumn])
                sedWavelength=float(bits[wavelengthColumn])
                sedFlux=float(bits[fluxColumn])
                
                row=self.ages.index(sedAge)
                column=self.wavelengths.index(sedWavelength)
                
                self.fluxGrid[row][column]=sedFlux

    def getSED(self, ageGyr, z = 0.0, normalise = False, label = None):
        """Extract a SED for given age. Do linear interpolation between models if necessary.
        
        @type ageGyr: float
        @param ageGyr: age of the SED in Gyr
        @type z: float
        @param z: redshift the SED from z = 0 to z = z
        @type normalise: bool
        @param normalise: normalise the SED to have area 1
        @rtype: L{SED} object
        @return: SED
        
        """
        
        if ageGyr in self.ages:
            
            flux=self.fluxGrid[self.ages.index(ageGyr)]
            sed=SED(self.wavelengths, flux, z = z, normalise = normalise, label = label)
            return sed
        
        else:
            
            # Use interpolation, iterating over each wavelength column
            flux=[]
            for i in range(len(self.wavelengths)):
                interpolator=interpolate.interp1d(self.ages, self.fluxGrid[:,i], kind='linear')
                sedFlux=interpolator(ageGyr)
                flux.append(sedFlux)
            sed=SED(self.wavelengths, flux, z = z, normalise = normalise, label = label)
            return sed

    def getColourEvolution(self, passband1, passband2, zFormation, zStepSize = 0.05, magType = "Vega"):
        """Calculates the evolution of the colour observed through passband1 - passband2 for the
        StellarPopulation with redshift, from z = 0 to z = zFormation.
        
        @type passband1: L{Passband} object
        @param passband1: filter passband through which to calculate the first magnitude
        @type passband2: L{Passband} object
        @param passband2: filter passband through which to calculate the second magnitude
        @type zFormation: float
        @param zFormation: formation redshift of the StellarPopulation
        @type zStepSize: float
        @param zStepSize: size of interval in z at which to calculate model colours
        @type magType: string
        @param magType: either "Vega" or "AB"
        @rtype: dictionary
        @return: dictionary of numpy.arrays in format {'z', 'colour'}
        
        """
       
        zSteps=int(math.ceil(zFormation/zStepSize))
        zData=[]
        colourData=[]
        for i in range(1, zSteps):
            zc=i*zStepSize
            age=astCalc.tl(zFormation)-astCalc.tl(zc)
            sed=self.getSED(age, z = zc)
            colour=sed.calcColour(passband1, passband2, magType = magType)
            zData.append(zc)
            colourData.append(colour)

        zData=numpy.array(zData)
        colourData=numpy.array(colourData)
        
        return {'z': zData, 'colour': colourData}
        
    def getMagEvolution(self, passband, magNormalisation, zNormalisation, zFormation, zStepSize = 0.05, 
                            onePlusZSteps = False, magType = "Vega"):
        """Calculates the evolution with redshift (from z = 0 to z = zFormation) of apparent magnitude
        in the observed frame through the passband for the StellarPopulation, normalised to magNormalisation 
        (apparent) at z = zNormalisation.
        
        @type passband: L{Passband} object
        @param passband: filter passband through which to calculate the magnitude
        @type magNormalisation: float
        @param magNormalisation: sets the apparent magnitude of the SED at zNormalisation
        @type zNormalisation: float
        @param zNormalisation: the redshift at which the magnitude normalisation is carried out
        @type zFormation: float
        @param zFormation: formation redshift of the StellarPopulation
        @type zStepSize: float
        @param zStepSize: size of interval in z at which to calculate model magnitudes
        @type onePlusZSteps: bool
        @param onePlusZSteps: if True, zSteps are (1+z)*zStepSize, otherwise zSteps are linear
        @type magType: string
        @param magType: either "Vega" or "AB"
        @rtype: dictionary
        @return: dictionary of numpy.arrays in format {'z', 'mag'}
        
        """
        
        # Count upwards in z steps as interpolation doesn't work if array ordered z decreasing
        zSteps=int(math.ceil(zFormation/zStepSize))
        zData=[]
        magData=[]
        absMagData=[]
        zc0=0.0
        for i in range(1, zSteps):
            if onePlusZSteps == False:
                zc=i*zStepSize
            else:
                zc=zc0+(1+zc0)*zStepSize
                zc0=zc
                if zc >= zFormation:
                    break
            age=astCalc.tl(zFormation)-astCalc.tl(zc)
            sed=self.getSED(age, z = zc)
            mag=sed.calcMag(passband, magType = magType, addDistanceModulus = True)
            zData.append(zc)
            magData.append(mag)
            absMagData.append(sed.calcMag(passband, addDistanceModulus = False))

        zData=numpy.array(zData)
        magData=numpy.array(magData)
        
        # Do the normalisation
        interpolator=interpolate.interp1d(zData, magData, kind='linear')
        modelNormMag=interpolator(zNormalisation)
        normConstant=magNormalisation-modelNormMag
        magData=magData+normConstant
        
        return {'z': zData, 'mag': magData}

    def calcEvolutionCorrection(self, zFrom, zTo, zFormation, passband, magType = "Vega"):
        """Calculates the evolution correction in magnitudes in the rest frame through the passband
        from redshift zFrom to redshift zTo, where the stellarPopulation is assumed to be formed
        at redshift zFormation.

        @type zFrom: float
        @param zFormation: redshift to evolution correct from
        @type zTo: float
        @param zTo: redshift to evolution correct to
        @type zFormation: float
        @param zFormation: formation redshift of the StellarPopulation
        @type passband: L{Passband} object
        @param passband: filter passband through which to calculate magnitude
        @type magType: string
        @param magType: either "Vega" or "AB"
        @rtype: float
        @return: evolution correction in magnitudes in the rest frame
        
        """
        
        ageFrom=astCalc.tl(zFormation)-astCalc.tl(zFrom)
        ageTo=astCalc.tl(zFormation)-astCalc.tl(zTo)
        
        fromSED=self.getSED(ageFrom)
        toSED=self.getSED(ageTo)
        
        fromMag=fromSED.calcMag(passband, magType = magType, addDistanceModulus = False)
        toMag=toSED.calcMag(passband, magType = magType, addDistanceModulus = False)
        
        return fromMag-toMag
        
#------------------------------------------------------------------------------------------------------------
class M05Model(StellarPopulation):
    """This class describes a Maraston 2005 stellar population model. To load a composite stellar population
    model (CSP) for a tau = 0.1 Gyr burst of star formation, solar metallicity, Salpeter IMF:
    
    m05csp = astSED.M05Model(M05_DIR+"/csp_e_0.10_z02_salp.sed_agb")
    
    where M05_DIR is set to point to the directory where the Maraston 2005 models are stored on your system.
    
    The file format of the Maraston 2005 simple stellar poulation (SSP) models is different to the file
    format used for the CSPs, and this needs to be specified using the fileType parameter. To load a SSP with
    solar metallicity, red horizontal branch morphology:
    
    m05ssp = astSED.M05Model(M05_DIR+"/sed.ssz002.rhb", fileType = "ssp")
    
    The wavelength units of SEDs from M05 models are Angstroms, with flux in units of erg/s/Angstrom.
    
    """
    def __init__(self, fileName, fileType = "csp"):
   
        self.modelFamily="M05"

        inFile=open(fileName, "r")
        lines=inFile.readlines()
        inFile.close()
        
        self.fileName=fileName

        if fileType == "csp":
            ageColumn=0
            wavelengthColumn=1
            fluxColumn=2
        elif fileType == "ssp":
            ageColumn=0
            wavelengthColumn=2
            fluxColumn=3
        else:
            raise Exception("fileType must be 'ssp' or 'csp'")
        
        # Extract a list of model ages and valid wavelengths from the file
        self.ages=[]
        self.wavelengths=[]
        for line in lines:
            if line[0] !="#" and len(line) > 3:
                bits=line.split()
                age=float(bits[ageColumn])
                wavelength=float(bits[wavelengthColumn])
                if age not in self.ages:
                    self.ages.append(age)
                if wavelength not in self.wavelengths:
                    self.wavelengths.append(wavelength)
        
        # Construct a grid of flux - rows correspond to each wavelength, columns to age
        self.fluxGrid=numpy.zeros([len(self.ages), len(self.wavelengths)])
        for line in lines:
            if line[0] !="#" and len(line) > 3:
                bits=line.split()
                sedAge=float(bits[ageColumn])
                sedWavelength=float(bits[wavelengthColumn])
                sedFlux=float(bits[fluxColumn])
                
                row=self.ages.index(sedAge)
                column=self.wavelengths.index(sedWavelength)
                
                self.fluxGrid[row][column]=sedFlux
    
#------------------------------------------------------------------------------------------------------------
class BC03Model(StellarPopulation):
    """This class describes a Bruzual & Charlot 2003 stellar population model, extracted from a GALAXEV .ised
    file using the galaxevpl program that is included in GALAXEV. The file format is white space delimited,
    with wavelength in the first column. Subsequent columns contain the model fluxes for SEDs of different
    ages, as specified when running galaxevpl. The age corresponding to each flux column is taken from the 
    comment line beginning "# Age (yr)", and is converted to Gyr.

    For example, to load a tau = 0.1 Gyr burst of star formation,  solar metallicity, Salpeter IMF model
    stored in a file (created by galaxevpl) called "csp_lr_solar_0p1Gyr.136":
    
    bc03model = BC03Model("csp_lr_solar_0p1Gyr.136")

    The wavelength units of SEDs from BC03 models are Angstroms. Flux is converted into units of 
    erg/s/Angstrom (the units in the files output by galaxevpl are LSun/Angstrom).

    """
    
    def __init__(self, fileName):
   
        self.modelFamily="BC03"
        self.fileName=fileName

        inFile=open(fileName, "r")
        lines=inFile.readlines()
        inFile.close()
        
        # Extract a list of model ages - BC03 ages are in years, so convert to Gyr
        self.ages=[]
        for line in lines:
            if line.find("# Age (yr)") != -1:
                rawAges=line[line.find("# Age (yr)")+10:].split()
                for age in rawAges:
                    self.ages.append(float(age)/1e9)
        
        # Extract a list of valid wavelengths from the file
        # If we have many ages in the file, this is more complicated...
        lambdaLinesCount=0
        startFluxDataLine=None
        for i in range(len(lines)):
            line=lines[i]
            if "# Lambda(A)" in line:
                lambdaLinesCount=lambdaLinesCount+1
            if line[0] != "#" and len(line) > 3 and startFluxDataLine == None:
                startFluxDataLine=i
        self.wavelengths=[]
        for i in range(startFluxDataLine, len(lines), lambdaLinesCount):
            line=lines[i]
            bits=line.split()
            self.wavelengths.append(float(bits[0]))        
        
        # Construct a grid of flux - rows correspond to each wavelength, columns to age
        self.fluxGrid=numpy.zeros([len(self.ages), len(self.wavelengths)])
        for i in range(startFluxDataLine, len(lines), lambdaLinesCount):
            line=lines[i]
            bits=[]
            for k in range(i, i+lambdaLinesCount):
                bits=bits+lines[k].split()           
            ageFluxes=bits[1:]
            sedWavelength=float(bits[0])
            column=self.wavelengths.index(sedWavelength)
            for row in range(len(ageFluxes)):
                sedFlux=float(ageFluxes[row])
                self.fluxGrid[row][column]=sedFlux

        # Convert flux into erg/s/Angstrom - native units of galaxevpl files are LSun/Angstrom
        self.fluxGrid=self.fluxGrid*3.826e33
        
#------------------------------------------------------------------------------------------------------------
class P09Model(StellarPopulation):
    """This class describes a Percival et al 2009 (BaSTI; http://albione.oa-teramo.inaf.it) stellar 
    population model. We assume that the synthetic spectra for each model are unpacked under the directory 
    pointed to by fileName.
    
    The wavelength units of SEDs from P09 models are converted to Angstroms. Flux is converted into units of 
    erg/s/Angstrom (the units in the BaSTI low-res spectra are 4.3607e-33 erg/s/m).
    
    """
    
    def __init__(self, fileName):
   
        self.modelFamily="P09"

        files=glob.glob(fileName+os.path.sep+"*.t??????")
        
        self.fileName=fileName

        # Map end of filenames to ages in Gyr
        extensionAgeMap={}
        self.ages=[]
        for f in files:
            ext=f.split(".")[-1]
            ageGyr=float(f[-5:])/1e3
            self.ages.append(ageGyr)
            extensionAgeMap[ext]=ageGyr
        self.ages.sort()
        
        # Construct a grid of flux - rows correspond to each wavelength, columns to age
        self.wavelengths=None
        self.fluxGrid=None
        for i in range(len(self.ages)):
            for e in extensionAgeMap.keys():
                if extensionAgeMap[e] == self.ages[i]:
                    inFileName=glob.glob(fileName+os.path.sep+"*."+e)[0]
                    inFile=open(inFileName, "r")
                    lines=inFile.readlines()
                    inFile.close()
                    wavelength=[]
                    flux=[]
                    for line in lines:
                        bits=line.split()
                        wavelength.append(float(bits[0])*10.0)  # units in file are nm, not angstroms
                        flux.append(float(bits[1]))
                    if self.wavelengths == None:
                        self.wavelengths=wavelength
                    if self.fluxGrid == None:
                        self.fluxGrid=numpy.zeros([len(self.ages), len(self.wavelengths)])
                    self.fluxGrid[i]=flux                    

        # Convert flux into erg/s/Angstrom - native units in BaSTI files are 4.3607e-33 erg/s/m
        self.fluxGrid=self.fluxGrid/4.3607e-33/1e10
        
#------------------------------------------------------------------------------------------------------------
def makeModelSEDDictList(modelList, z, passbandsList, labelsList = [], EBMinusVList = [0.0], forceYoungerThanUniverse = True):
    """This routine makes a list of SEDDict dictionaries (see L{mags2SEDDict}) for fitting using 
    L{fitSEDDict}. This speeds up the fitting as this allows us to calculate model SED magnitudes only once, 
    if all objects to be fitted are at the same redshift. We add some meta data to the modelSEDDicts (e.g.
    the model file names).
        
    The effect of extinction by dust (assuming the Calzetti et al. 2000 law) can be included by giving a list 
    of E(B-V) values.
    
    If forceYoungerThanUniverse == True, ages which are older than the universe at the given z will not be
    included.
    
    @type modelList: list of L{StellarPopulation} model objects
    @param modelList: list of StellarPopulation models to include
    @type z: float
    @param z: redshift to apply to all stellar population models in modelList
    @type EBMinusVList: list
    @param EBMinusVList: list of E(B-V) extinction values to apply to all models, in magnitudes
    @type labelsList: list
    @param labelsList: optional list used for labelling passbands in output SEDDicts
    @type forceYoungerThanUniverse: bool
    @param forceYoungerThanUniverse: if True, do not allow models that exceed the age of the universe at z
    @rtype: list
    @return: list of dictionaries containing model fluxes, to be used as input to L{fitSEDDict}.
    
    """
    
    # Otherwise if this is the case we won't actually make any model SEDDicts ...
    if EBMinusVList == []:
        EBMinusVList=[0.0]
        
    modelSEDDictList=[]
    for m in range(len(modelList)):
        testAges=numpy.array(modelList[m].ages)
        if forceYoungerThanUniverse == True:
            testAges=testAges[numpy.logical_and(numpy.less(testAges, astCalc.tz(z)), numpy.greater(testAges, 0))]
        for t in testAges:
            s=modelList[m].getSED(t, z = z, label=modelList[m].fileName+" - age="+str(t)+" Gyr")
            for EBMinusV in EBMinusVList:
                try:
                    s.extinctionCalzetti(EBMinusV)
                except:
                    raise Exception("Model %s has a wavelength range that doesn't cover ~1200-22000 Angstroms" % (modelList[m].fileName))
                modelSEDDict=s.getSEDDict(passbandsList)
                modelSEDDict['labels']=labelsList
                modelSEDDict['E(B-V)']=EBMinusV
                modelSEDDict['ageGyr']=t
                modelSEDDict['z']=z
                modelSEDDict['fileName']=modelList[m].fileName               
                modelSEDDict['modelListIndex']=m
                modelSEDDictList.append(modelSEDDict)
    
    return modelSEDDictList
    
#------------------------------------------------------------------------------------------------------------
def fitSEDDict(SEDDict, modelSEDDictList):
    """Fits the given SED dictionary (made using L{mags2SEDDict}) with the given list of model SED 
    dictionaries. The latter should be made using L{makeModelSEDDictList}, and entries for fluxes should
    correspond directly between the model and SEDDict.
           
    Returns a dictionary with best fit values.
    
    @type SEDDict: dictionary, in format of L{mags2SEDDict}
    @param SEDDict: dictionary of observed fluxes and uncertainties, in format of L{mags2SEDDict}
    @type modelSEDDictList: list of dictionaries, in format of L{makeModelSEDDictList}
    @param modelSEDDictList: list of dictionaries containing fluxes of models to be fitted to the observed
    fluxes listed in the SEDDict. This should be made using L{makeModelSEDDictList}.
    @rtype: dictionary
    @return: results of the fitting - keys: 
             - 'minChiSq': minimum chi squared value of best fit
             - 'chiSqContrib': corresponding contribution at each passband to the minimum chi squared value
             - 'ageGyr': the age in Gyr of the best fitting model
             - 'modelFileName': the file name of the stellar population model corresponding to the best fit
             - 'modelListIndex': the index of the best fitting model in the input modelSEDDictList
             - 'norm': the normalisation that the best fit model should be multiplied by to match the SEDDict
             - 'z': the redshift of the best fit model
             - 'E(B-V)': the extinction, E(B-V), in magnitudes, of the best fit model
    
    """
    
    modelFlux=[]
    for modelSEDDict in modelSEDDictList:
        modelFlux.append(modelSEDDict['flux'])
    modelFlux=numpy.array(modelFlux)    
    sedFlux=numpy.array([SEDDict['flux']]*len(modelSEDDictList))
    sedFluxErr=numpy.array([SEDDict['fluxErr']]*len(modelSEDDictList))

    # Analytic expression below is for normalisation at minimum chi squared (see note book)
    norm=numpy.sum((modelFlux*sedFlux)/(sedFluxErr**2), axis=1)/numpy.sum(modelFlux**2/sedFluxErr**2, axis=1)
    norms=numpy.array([norm]*modelFlux.shape[1]).transpose()
    chiSq=numpy.sum(((sedFlux-norms*modelFlux)**2)/sedFluxErr**2, axis=1)
    chiSq[numpy.isnan(chiSq)]=1e6   # throw these out, should check this out and handle more gracefully
    minChiSq=chiSq.min()
    bestMatchIndex=numpy.equal(chiSq, minChiSq).nonzero()[0][0]
    bestNorm=norm[bestMatchIndex]
    bestChiSq=minChiSq
    bestChiSqContrib=((sedFlux[bestMatchIndex]-norms[bestMatchIndex]*modelFlux[bestMatchIndex])**2)\
                        /sedFluxErr[bestMatchIndex]**2
    
    resultsDict={'minChiSq': bestChiSq, 
                 'chiSqContrib': bestChiSqContrib,
                 'allChiSqValues': chiSq,
                 'ageGyr': modelSEDDictList[bestMatchIndex]['ageGyr'], 
                 'modelFileName': modelSEDDictList[bestMatchIndex]['fileName'],
                 'modelListIndex': modelSEDDictList[bestMatchIndex]['modelListIndex'],
                 'norm': bestNorm, 
                 'z': modelSEDDictList[bestMatchIndex]['z'], 
                 'E(B-V)': modelSEDDictList[bestMatchIndex]['E(B-V)']}
    
    return resultsDict
    
#------------------------------------------------------------------------------------------------------------
def mags2SEDDict(ABMags, ABMagErrs, passbands):
    """Takes a set of corresponding AB magnitudes, uncertainties, and passbands, and
    returns a dictionary with keys 'flux', 'fluxErr', 'wavelength'. Fluxes are in units of 
    erg/s/cm^2/Angstrom, wavelength in Angstroms. These dictionaries are the staple diet of the
    L{fitSEDDict} routine.
    
    @type ABMags: list or numpy array
    @param ABMags: AB magnitudes, specified in corresponding order to passbands and ABMagErrs
    @type ABMagErrs: list or numpy array
    @param ABMagErrs: AB magnitude errors, specified in corresponding order to passbands and ABMags
    @type passbands: list of L{Passband} objects
    @param passbands: passband objects, specified in corresponding order to ABMags and ABMagErrs
    @rtype: dictionary
    @return: dictionary with keys {'flux', 'fluxErr', 'wavelength'}, suitable for input to L{fitSEDDict}

    """
    
    flux=[]
    fluxErr=[]
    wavelength=[]
    for m, e, p in zip(ABMags, ABMagErrs, passbands):
        f, err=mag2Flux(m, e, p)
        flux.append(f)
        fluxErr.append(err)
        wavelength.append(p.effectiveWavelength())
        
    SEDDict={}
    SEDDict['flux']=numpy.array(flux)
    SEDDict['fluxErr']=numpy.array(fluxErr)
    SEDDict['wavelength']=numpy.array(wavelength)
    
    return SEDDict
    
#------------------------------------------------------------------------------------------------------------
def mag2Flux(ABMag, ABMagErr, passband):
    """Converts given AB magnitude and uncertainty into flux, in erg/s/cm^2/Angstrom.
    
    @type ABMag: float
    @param ABMag: magnitude on AB system in passband
    @type ABMagErr: float
    @param ABMagErr: uncertainty in AB magnitude in passband
    @type passband: L{Passband} object
    @param passband: L{Passband} object at which ABMag was measured
    @rtype: list
    @return: [flux, fluxError], in units of erg/s/cm^2/Angstrom
    
    """
    
    fluxJy=(10**23.0)*10**(-(ABMag+48.6)/2.5)   # AB mag
    aLambda=3e-13 # for conversion to erg s-1 cm-2 angstrom-1 with lambda in microns
    effLMicron=passband.effectiveWavelength()*(1e-10/1e-6)
    fluxWLUnits=aLambda*fluxJy/effLMicron**2
        
    fluxJyErr=(10**23.0)*10**(-(ABMag-ABMagErr+48.6)/2.5)   # AB mag
    fluxWLUnitsErr=aLambda*fluxJyErr/effLMicron**2
    fluxWLUnitsErr=fluxWLUnitsErr-fluxWLUnits

    return [fluxWLUnits, fluxWLUnitsErr]

#------------------------------------------------------------------------------------------------------------
def flux2Mag(flux, fluxErr, passband):
    """Converts given flux and uncertainty in erg/s/cm^2/Angstrom into AB magnitudes.
    
    @type flux: float
    @param flux: flux in erg/s/cm^2/Angstrom in passband
    @type fluxErr: float
    @param fluxErr: uncertainty in flux in passband, in erg/s/cm^2/Angstrom
    @type passband: L{Passband} object
    @param passband: L{Passband} object at which ABMag was measured
    @rtype: list
    @return: [ABMag, ABMagError], in AB magnitudes
    
    """

    # aLambda = 3x10-5 for effective wavelength in angstroms
    aLambda=3e-13 # for conversion to erg s-1 cm-2 angstrom-1 with lambda in microns
    effLMicron=passband.effectiveWavelength()*(1e-10/1e-6)

    fluxJy=(flux*effLMicron**2)/aLambda
    mag=-2.5*numpy.log10(fluxJy/10**23)-48.6
    
    fluxErrJy=(fluxErr*effLMicron**2)/aLambda
    magErr=mag-(-2.5*numpy.log10((fluxJy+fluxErrJy)/10**23)-48.6)
    
    return [mag, magErr]

#------------------------------------------------------------------------------------------------------------
def mag2Jy(ABMag):
    """Converts an AB magnitude into flux density in Jy
    
    @type ABMag: float
    @param ABMag: AB magnitude
    @rtype: float
    @return: flux density in Jy
    
    """
    
    fluxJy=((10**23)*10**(-(float(ABMag)+48.6)/2.5))
    
    return fluxJy


#------------------------------------------------------------------------------------------------------------
def Jy2Mag(fluxJy):
    """Converts flux density in Jy into AB magnitude
    
    @type fluxJy: float
    @param fluxJy: flux density in Jy
    @rtype: float
    @return: AB magnitude
    
    """
        
    ABMag=-2.5*(numpy.log10(fluxJy)-23.0)-48.6
    
    return ABMag
    
#------------------------------------------------------------------------------------------------------------
# Data
VEGA=VegaSED()

# AB SED has constant flux density 3631 Jy
AB=SED(wavelength = numpy.logspace(1, 8, 1e5), flux = numpy.ones(1000000))
AB.flux=(3e-5*3631)/(AB.wavelength**2)
AB.z0flux=AB.flux[:]

# Solar SED from HST CALSPEC (http://www.stsci.edu/hst/observatory/cdbs/calspec.html)
SOL=SED()
SOL.loadFromFile(astLib.__path__[0]+os.path.sep+"data"+os.path.sep+"sun_reference_stis_001.ascii")

