.. image:: https://readthedocs.org/projects/astlib/badge/?version=latest

**astLib** provides some tools for research astronomers who use Python.

* **Documentation:** https://astlib.readthedocs.io
* **License:** `LGPL <LICENSE>`_
* **Authors:** Matt Hilton & Steven Boada
* **Installation:** ``pip install astLib``
* **Support:** Please use the `GitHub issues page <https://github.com/mattyowl/astLib/issues>`_,
  and/or contact `Matt Hilton <mailto:matt.hilton@mykolab.com>`_.

**astLib** is divided into several modules:

* astCalc   (general calculations, e.g. luminosity distance etc.)
* astCoords (coordinate conversions etc.)
* astImages (clip sections from .fits etc.) 
* astPlots  (provides a flexible image plot class, e.g. plot image with catalogue objects overlaid)
* astSED    (calculate colours, magnitudes from stellar population models or spectral templates, fit photometric observations using stellar population models etc.)
* astStats  (statistics, e.g. biweight location/scale estimators etc.)
* astWCS    (routines for using FITS World Coordinate System information)

The astWCS module is a higher level interface to PyWCSTools, a simple SWIG (http://www.swig.org) wrapping 
of some of the routines from WCSTools by Jessica Mink (http://tdc-www.harvard.edu/software/wcstools/). It is 
used by some routines in astCoords, astImages and astPlots.

The goal of **astLib** was to provide features useful to astronomers that are not included in the scipy
(http://scipy.org), numpy (http://numpy.scipy.org) or matplotlib (http://matplotlib.sourceforge.net) modules 
on which astLib depends. For a far more extensive set of Python astronomy modules, see astropy 
(http://www.astropy.org/).

Some scripts using **astLib** can be found in the examples/ folder provided with the source code distribution.
