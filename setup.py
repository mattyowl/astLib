import os
import glob
import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
import sysconfig
from pkg_resources import require

topDir = os.getcwd()
sourceDir = "PyWCSTools"+os.path.sep+"wcssubs-3.9.5"+os.path.sep

#oFiles=glob.glob(sourceDir+"*.o")
#print oFiles
oFiles = ['PyWCSTools/wcssubs-3.9.5/cel.o', 'PyWCSTools/wcssubs-3.9.5/wcs.o',
    'PyWCSTools/wcssubs-3.9.5/proj.o', 'PyWCSTools/wcssubs-3.9.5/distort.o',
    'PyWCSTools/wcssubs-3.9.5/wcsinit.o', 'PyWCSTools/wcssubs-3.9.5/wcslib.o',
    'PyWCSTools/wcssubs-3.9.5/poly.o', 'PyWCSTools/wcssubs-3.9.5/platepos.o',
    'PyWCSTools/wcssubs-3.9.5/zpxpos.o', 'PyWCSTools/wcssubs-3.9.5/iget.o',
    'PyWCSTools/wcssubs-3.9.5/imio.o', 'PyWCSTools/wcssubs-3.9.5/dsspos.o',
    'PyWCSTools/wcssubs-3.9.5/tnxpos.o', 'PyWCSTools/wcssubs-3.9.5/wcscon.o',
    'PyWCSTools/wcssubs-3.9.5/fitsfile.o',
    'PyWCSTools/wcssubs-3.9.5/dateutil.o',
    'PyWCSTools/wcssubs-3.9.5/imhfile.o', 'PyWCSTools/wcssubs-3.9.5/lin.o',
    'PyWCSTools/wcssubs-3.9.5/fileutil.o',
    'PyWCSTools/wcssubs-3.9.5/wcstrig.o',
    'PyWCSTools/wcssubs-3.9.5/sph.o',
    'PyWCSTools/wcssubs-3.9.5/worldpos.o', 'PyWCSTools/wcssubs-3.9.5/hget.o',
    'PyWCSTools/wcssubs-3.9.5/hput.o']

exampleScripts = glob.glob("scripts"+os.path.sep+"*.py")

class build_PyWCSTools_ext(build_ext):

    def build_extensions(self):

        os.chdir(sourceDir)
        cc=setuptools._distutils.ccompiler.new_compiler(setuptools._distutils.ccompiler.get_default_compiler())
        cc.compiler_so=[sysconfig.get_config_var('CC')]+sysconfig.get_config_var('CFLAGS').split()+sysconfig.get_config_var('CFLAGSFORSHARED').split()

        # Suppress warnings from compiling WCSTools wcssubs-3.9.5
        if "-Wstrict-prototypes" in cc.compiler_so:
            cc.compiler_so.pop(cc.compiler_so.index("-Wstrict-prototypes"))
        if "-Wall" in cc.compiler_so:
            cc.compiler_so.pop(cc.compiler_so.index("-Wall"))
        # For recent macOS
        if "-Wno-error=implicit-function-declaration" in cc.compiler_so:
            cc.compiler_so.pop(cc.compiler_so.index("-Wno-error=implicit-function-declaration"))

        WCSToolsCFiles = glob.glob("*.c")
        WCSToolsCFiles.pop(WCSToolsCFiles.index("wcs_wrap.c"))
        WCSToolsCFiles.pop(WCSToolsCFiles.index("wcscon_wrap.c"))
        cc.compile(WCSToolsCFiles)
        os.chdir(topDir)

        build_ext.build_extensions(self)

setup(name='astLib',
    version='0.11.9',
    url='https://astlib.readthedocs.io',
    author='Matt Hilton & Steven Boada',
    author_email='matt.hilton@mykolab.com',
    classifiers=['Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
            'Natural Language :: English',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Astronomy',
            'Topic :: Software Development :: Libraries'],
    description='A set of python modules for producing simple plots, statistics, common calculations, coordinate conversions, and manipulating FITS images with World Coordinate System (WCS) information.',
    long_description="""astLib is a set of Python modules that provides some tools for research astronomers. It can be
used for simple plots, statistics, common calculations, coordinate conversions, and manipulating FITS images
with World Coordinate System (WCS) information through PyWCSTools - a simple wrapping of WCSTools by Jessica Mink.
PyWCSTools is distributed (and developed) as part of astLib.""",
    packages=['astLib', 'PyWCSTools'],
    package_data={'astLib': ['data/*']},
    cmdclass={"build_ext": build_PyWCSTools_ext},
    scripts=exampleScripts,
    ext_modules=[
        Extension('PyWCSTools._wcscon', [sourceDir+"wcscon_wrap.c"],
        extra_objects=oFiles),
        Extension('PyWCSTools._wcs', [sourceDir+"wcs_wrap.c"],
        extra_objects=oFiles)
    ]
	)
