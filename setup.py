import os
import glob
import setuptools
from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.extension import Extension
import sysconfig
from pkg_resources import require

topDir = os.getcwd()
sourceDir = "PyWCSTools"+os.path.sep+"wcssubs-3.9.7"+os.path.sep

#oFiles=glob.glob(sourceDir+"*.o")
#print oFiles
oFiles = ['PyWCSTools/wcssubs-3.9.7/cel.o', 'PyWCSTools/wcssubs-3.9.7/wcs.o',
    'PyWCSTools/wcssubs-3.9.7/proj.o', 'PyWCSTools/wcssubs-3.9.7/distort.o',
    'PyWCSTools/wcssubs-3.9.7/wcsinit.o', 'PyWCSTools/wcssubs-3.9.7/wcslib.o',
    'PyWCSTools/wcssubs-3.9.7/poly.o', 'PyWCSTools/wcssubs-3.9.7/platepos.o',
    'PyWCSTools/wcssubs-3.9.7/zpxpos.o', 'PyWCSTools/wcssubs-3.9.7/iget.o',
    'PyWCSTools/wcssubs-3.9.7/imio.o', 'PyWCSTools/wcssubs-3.9.7/dsspos.o',
    'PyWCSTools/wcssubs-3.9.7/tnxpos.o', 'PyWCSTools/wcssubs-3.9.7/wcscon.o',
    'PyWCSTools/wcssubs-3.9.7/fitsfile.o',
    'PyWCSTools/wcssubs-3.9.7/dateutil.o',
    'PyWCSTools/wcssubs-3.9.7/imhfile.o', 'PyWCSTools/wcssubs-3.9.7/lin.o',
    'PyWCSTools/wcssubs-3.9.7/fileutil.o',
    'PyWCSTools/wcssubs-3.9.7/wcstrig.o',
    'PyWCSTools/wcssubs-3.9.7/sph.o',
    'PyWCSTools/wcssubs-3.9.7/worldpos.o', 'PyWCSTools/wcssubs-3.9.7/hget.o',
    'PyWCSTools/wcssubs-3.9.7/hput.o']

exampleScripts = glob.glob("scripts"+os.path.sep+"*.py")

class build_PyWCSTools_ext(build_ext):

    def build_extensions(self):

        os.chdir(sourceDir)
        cc=setuptools._distutils.ccompiler.new_compiler(setuptools._distutils.ccompiler.get_default_compiler())
        cc.compiler_so=sysconfig.get_config_var('CC').split()+sysconfig.get_config_var('CFLAGS').split()+sysconfig.get_config_var('CFLAGSFORSHARED').split()
        if '-fPIC' not in cc.compiler_so:
            cc.compiler_so.append('-fPIC') # Needed by Ilifu but not picked automatically

        # Suppress warnings from compiling WCSTools wcssubs-3.9.7
        if "-Wstrict-prototypes" in cc.compiler_so:
            cc.compiler_so.pop(cc.compiler_so.index("-Wstrict-prototypes"))
        if "-Wall" in cc.compiler_so:
            cc.compiler_so.pop(cc.compiler_so.index("-Wall"))
        # For recent macOS
        # if "-Wno-error=implicit-function-declaration" in cc.compiler_so:
        #     cc.compiler_so.pop(cc.compiler_so.index("-Wno-error=implicit-function-declaration"))

        WCSToolsCFiles = glob.glob("*.c")
        WCSToolsCFiles.pop(WCSToolsCFiles.index("wcs_wrap.c"))
        WCSToolsCFiles.pop(WCSToolsCFiles.index("wcscon_wrap.c"))
        cc.compile(WCSToolsCFiles)
        os.chdir(topDir)

        build_ext.build_extensions(self)

setup(name='astLib',
    version='0.13.1',
    packages=['astLib', 'PyWCSTools'],
    package_data={'astLib': ['data/*']},
    cmdclass={"build_ext": build_PyWCSTools_ext},
    scripts=exampleScripts,
    ext_modules=[
        # Extension('PyWCSTools._wcscon', [sourceDir+"wcscon_wrap.c"],
        # extra_objects=oFiles),
        Extension('PyWCSTools._wcs', [sourceDir+"wcs_wrap.c"],
        extra_objects=oFiles)
    ]
	)
