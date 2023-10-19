.. _UsagePage:
  
=============================
Usage
=============================

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


Known issues (may no longer be relevant)
========================================

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
