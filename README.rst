.. image:: https://readthedocs.org/projects/astlib/badge/?version=latest

astLib provides some tools for research astronomers who use Python. It is divided into several modules:

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

The goal of astLib was to provide features useful to astronomers that are not included in the scipy 
(http://scipy.org), numpy (http://numpy.scipy.org) or matplotlib (http://matplotlib.sourceforge.net) modules 
on which astLib depends. For a far more extensive set of Python astronomy modules, see astropy 
(http://www.astropy.org/).

Some scripts using astLib can be found in the examples/ folder provided with the source code distribution.


Software needed
===============

astLib requires:

* Python (tested on versions 3.6+)
* Astropy - http://www.astropy.org (tested on version 3.2.1)
* Numpy - http://numpy.scipy.org (tested on version 1.18.1)
* SciPy - http://scipy.org (tested on version 1.3.1)
* Matplotlib - http://matplotlib.sourceforge.net (tested on version 3.1.1)

Optional:
   
* Python Imaging Library - http://www.pythonware.com/products/pil (tested on version 1.1.7)

Other versions of the software listed above are likely to work.


Installation
============

You can install astLib via pip:

.. code-block::

   pip install astLib --user


You may also install using the standard ``setup.py`` script, e.g., as root:

.. code-block::

   sudo python setup.py install


Alternatively, 

.. code-block::

   python setup.py install --user


will install ``astLib`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

You can also use the ``--prefix`` option, e.g.,

.. code-block::

   python setup.py install --prefix=$HOME/local


and then add, e.g., ``$HOME/local/lib/python3.6/site-packages`` to 
$PYTHONPATH (adjust the path according to your Python version number).

.. code-block::

   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH


Installation on recent versions of macOS
========================================

Some users have reported that the standard method for installing ``astLib`` does not work on recent versions
of macOS (e.g., Big Sur), due to the default compiler flags. The current workaround for this is to install
using:
  
.. code-block::

   CFLAGS="-Wno-error=implicit-function-declaration" python setup.py install
   

Thanks to Michael Cowley and Stefano Covino for helping to resolve this issue.


Usage
=====

To access the routines in the astLib modules, simply:

.. code-block::

   from astLib import astCalc
   from astLib import astCoords
   from astLib import astWCS


etc.

The astWCS module currently provides access to what are (I think) the most commonly needed WCS information 
and functions (such as converting between pixel and WCS coordinates etc.). However, should you wish to 
access the wrapped WCSTools routines themselves directly: 

.. code-block::

   from PyWCSTools import wcs
   from PyWCSTools import wcscon

etc.

Note that PyWCSTools only includes some functions from wcs.c and wcscon.c at present. For examples of usage, 
look at the Python code for the astLib.astWCS module. Documentation for the WCSTools routines can be found 
here: http://tdc-www.harvard.edu/software/wcstools/subroutines/libwcs.wcs.html.

As of version 0.11.x+, by default the ``astWCS.WCS`` class is using the ``astropy.wcs`` module instead of
PyWCSTools (this allows one to benefit from some features of ``astropy.wcs`` without having to re-write
code based on ``astWCS.WCS``). To use PyWCSTools instead, set ``useAstropyWCS = False`` when creating a
``WCS`` object.


Known issues
============

This may no longer apply, but just in case...

Recent versions of matplotlib (on which astLib depends) now use locale information. On systems where the
decimal point separator is not '.' (e.g. Germany), the astWCS coordinate conversions routines will give
strange results if this is not accounted for. As of version 0.3.0, the astWCS module will detect if this is 
the case and print a warning message to the console.

The workaround for this issue is to add the following after importing any python modules that expicitly set 
the locale (such as matplotlib):

.. code-block::
    
    import locale
    locale.setlocale(locale.LC_NUMERIC, 'C')"

Thanks to Markus Demleitner for pointing this out.


Documentation
=============

Documentation is available on the web at http://astlib.readthedocs.io.


Bugs
====

Please email bug reports to matt.hilton@mykolab.com, and/or use the `GitHub issues page <https://github.com/mattyowl/astLib/issues>`_. 

Please include details of your operating system, python version, and versions of the python packages
required by astLib that you have installed on your machine. For any WCS-related bugs, it would be helpful 
if you could also include the image header as a text file so that I can reproduce them easily. 
