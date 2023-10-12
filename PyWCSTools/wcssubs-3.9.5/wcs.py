# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _wcs
else:
    import _wcs

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)



def new_doubleArray(nelements):
    r"""new_doubleArray(nelements) -> double *"""
    return _wcs.new_doubleArray(nelements)

def delete_doubleArray(ary):
    r"""delete_doubleArray(ary)"""
    return _wcs.delete_doubleArray(ary)

def doubleArray_getitem(ary, index):
    r"""doubleArray_getitem(ary, index) -> double"""
    return _wcs.doubleArray_getitem(ary, index)

def doubleArray_setitem(ary, index, value):
    r"""doubleArray_setitem(ary, index, value)"""
    return _wcs.doubleArray_setitem(ary, index, value)

def wcsinit(hstring):
    r"""wcsinit(hstring) -> WorldCoor"""
    return _wcs.wcsinit(hstring)

def wcsxinit(cra, cdec, secpix, xrpix, yrpix, nxpix, nypix, rotate, equinox, epoch, proj):
    r"""wcsxinit(cra, cdec, secpix, xrpix, yrpix, nxpix, nypix, rotate, equinox, epoch, proj) -> WorldCoor"""
    return _wcs.wcsxinit(cra, cdec, secpix, xrpix, yrpix, nxpix, nypix, rotate, equinox, epoch, proj)

def wcskinit(nxpix, nypix, ctype1, ctype2, crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2, crota, equinox, epoch):
    r"""wcskinit(nxpix, nypix, ctype1, ctype2, crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2, crota, equinox, epoch) -> WorldCoor"""
    return _wcs.wcskinit(nxpix, nypix, ctype1, ctype2, crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2, crota, equinox, epoch)

def iswcs(wcs):
    r"""iswcs(wcs) -> int"""
    return _wcs.iswcs(wcs)

def nowcs(wcs):
    r"""nowcs(wcs) -> int"""
    return _wcs.nowcs(wcs)

def wcs2pix(wcs, xpos, ypos):
    r"""wcs2pix(wcs, xpos, ypos)"""
    return _wcs.wcs2pix(wcs, xpos, ypos)

def pix2wcs(wcs, xpix, ypix):
    r"""pix2wcs(wcs, xpix, ypix)"""
    return _wcs.pix2wcs(wcs, xpix, ypix)

def wcscent(wcs):
    r"""wcscent(wcs)"""
    return _wcs.wcscent(wcs)

def getradecsys(wcs):
    r"""getradecsys(wcs) -> char *"""
    return _wcs.getradecsys(wcs)

def wcsoutinit(wcs, coorsys):
    r"""wcsoutinit(wcs, coorsys)"""
    return _wcs.wcsoutinit(wcs, coorsys)

def wcsininit(wcs, coorsys):
    r"""wcsininit(wcs, coorsys)"""
    return _wcs.wcsininit(wcs, coorsys)

def getwcsout(wcs):
    r"""getwcsout(wcs) -> char *"""
    return _wcs.getwcsout(wcs)

def getwcsin(wcs):
    r"""getwcsin(wcs) -> char *"""
    return _wcs.getwcsin(wcs)

def wcssize(wcs):
    r"""wcssize(wcs)"""
    return _wcs.wcssize(wcs)

def wcsfull(wcs):
    r"""wcsfull(wcs)"""
    return _wcs.wcsfull(wcs)
class WorldCoor(object):
    r"""Proxy of C WorldCoor struct."""

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    xref = property(_wcs.WorldCoor_xref_get, _wcs.WorldCoor_xref_set, doc=r"""xref""")
    yref = property(_wcs.WorldCoor_yref_get, _wcs.WorldCoor_yref_set, doc=r"""yref""")
    xrefpix = property(_wcs.WorldCoor_xrefpix_get, _wcs.WorldCoor_xrefpix_set, doc=r"""xrefpix""")
    yrefpix = property(_wcs.WorldCoor_yrefpix_get, _wcs.WorldCoor_yrefpix_set, doc=r"""yrefpix""")
    xinc = property(_wcs.WorldCoor_xinc_get, _wcs.WorldCoor_xinc_set, doc=r"""xinc""")
    yinc = property(_wcs.WorldCoor_yinc_get, _wcs.WorldCoor_yinc_set, doc=r"""yinc""")
    rot = property(_wcs.WorldCoor_rot_get, _wcs.WorldCoor_rot_set, doc=r"""rot""")
    cd = property(_wcs.WorldCoor_cd_get, _wcs.WorldCoor_cd_set, doc=r"""cd""")
    dc = property(_wcs.WorldCoor_dc_get, _wcs.WorldCoor_dc_set, doc=r"""dc""")
    equinox = property(_wcs.WorldCoor_equinox_get, _wcs.WorldCoor_equinox_set, doc=r"""equinox""")
    epoch = property(_wcs.WorldCoor_epoch_get, _wcs.WorldCoor_epoch_set, doc=r"""epoch""")
    nxpix = property(_wcs.WorldCoor_nxpix_get, _wcs.WorldCoor_nxpix_set, doc=r"""nxpix""")
    nypix = property(_wcs.WorldCoor_nypix_get, _wcs.WorldCoor_nypix_set, doc=r"""nypix""")
    plate_ra = property(_wcs.WorldCoor_plate_ra_get, _wcs.WorldCoor_plate_ra_set, doc=r"""plate_ra""")
    plate_dec = property(_wcs.WorldCoor_plate_dec_get, _wcs.WorldCoor_plate_dec_set, doc=r"""plate_dec""")
    plate_scale = property(_wcs.WorldCoor_plate_scale_get, _wcs.WorldCoor_plate_scale_set, doc=r"""plate_scale""")
    x_pixel_offset = property(_wcs.WorldCoor_x_pixel_offset_get, _wcs.WorldCoor_x_pixel_offset_set, doc=r"""x_pixel_offset""")
    y_pixel_offset = property(_wcs.WorldCoor_y_pixel_offset_get, _wcs.WorldCoor_y_pixel_offset_set, doc=r"""y_pixel_offset""")
    x_pixel_size = property(_wcs.WorldCoor_x_pixel_size_get, _wcs.WorldCoor_x_pixel_size_set, doc=r"""x_pixel_size""")
    y_pixel_size = property(_wcs.WorldCoor_y_pixel_size_get, _wcs.WorldCoor_y_pixel_size_set, doc=r"""y_pixel_size""")
    ppo_coeff = property(_wcs.WorldCoor_ppo_coeff_get, _wcs.WorldCoor_ppo_coeff_set, doc=r"""ppo_coeff""")
    x_coeff = property(_wcs.WorldCoor_x_coeff_get, _wcs.WorldCoor_x_coeff_set, doc=r"""x_coeff""")
    y_coeff = property(_wcs.WorldCoor_y_coeff_get, _wcs.WorldCoor_y_coeff_set, doc=r"""y_coeff""")
    xpix = property(_wcs.WorldCoor_xpix_get, _wcs.WorldCoor_xpix_set, doc=r"""xpix""")
    ypix = property(_wcs.WorldCoor_ypix_get, _wcs.WorldCoor_ypix_set, doc=r"""ypix""")
    zpix = property(_wcs.WorldCoor_zpix_get, _wcs.WorldCoor_zpix_set, doc=r"""zpix""")
    xpos = property(_wcs.WorldCoor_xpos_get, _wcs.WorldCoor_xpos_set, doc=r"""xpos""")
    ypos = property(_wcs.WorldCoor_ypos_get, _wcs.WorldCoor_ypos_set, doc=r"""ypos""")
    crpix = property(_wcs.WorldCoor_crpix_get, _wcs.WorldCoor_crpix_set, doc=r"""crpix""")
    crval = property(_wcs.WorldCoor_crval_get, _wcs.WorldCoor_crval_set, doc=r"""crval""")
    cdelt = property(_wcs.WorldCoor_cdelt_get, _wcs.WorldCoor_cdelt_set, doc=r"""cdelt""")
    pc = property(_wcs.WorldCoor_pc_get, _wcs.WorldCoor_pc_set, doc=r"""pc""")
    projp = property(_wcs.WorldCoor_projp_get, _wcs.WorldCoor_projp_set, doc=r"""projp""")
    pvfail = property(_wcs.WorldCoor_pvfail_get, _wcs.WorldCoor_pvfail_set, doc=r"""pvfail""")
    projppv = property(_wcs.WorldCoor_projppv_get, _wcs.WorldCoor_projppv_set, doc=r"""projppv""")
    inv_x = property(_wcs.WorldCoor_inv_x_get, _wcs.WorldCoor_inv_x_set, doc=r"""inv_x""")
    inv_y = property(_wcs.WorldCoor_inv_y_get, _wcs.WorldCoor_inv_y_set, doc=r"""inv_y""")
    longpole = property(_wcs.WorldCoor_longpole_get, _wcs.WorldCoor_longpole_set, doc=r"""longpole""")
    latpole = property(_wcs.WorldCoor_latpole_get, _wcs.WorldCoor_latpole_set, doc=r"""latpole""")
    rodeg = property(_wcs.WorldCoor_rodeg_get, _wcs.WorldCoor_rodeg_set, doc=r"""rodeg""")
    imrot = property(_wcs.WorldCoor_imrot_get, _wcs.WorldCoor_imrot_set, doc=r"""imrot""")
    pa_north = property(_wcs.WorldCoor_pa_north_get, _wcs.WorldCoor_pa_north_set, doc=r"""pa_north""")
    pa_east = property(_wcs.WorldCoor_pa_east_get, _wcs.WorldCoor_pa_east_set, doc=r"""pa_east""")
    radvel = property(_wcs.WorldCoor_radvel_get, _wcs.WorldCoor_radvel_set, doc=r"""radvel""")
    zvel = property(_wcs.WorldCoor_zvel_get, _wcs.WorldCoor_zvel_set, doc=r"""zvel""")
    zpzd = property(_wcs.WorldCoor_zpzd_get, _wcs.WorldCoor_zpzd_set, doc=r"""zpzd""")
    zpr = property(_wcs.WorldCoor_zpr_get, _wcs.WorldCoor_zpr_set, doc=r"""zpr""")
    imflip = property(_wcs.WorldCoor_imflip_get, _wcs.WorldCoor_imflip_set, doc=r"""imflip""")
    prjcode = property(_wcs.WorldCoor_prjcode_get, _wcs.WorldCoor_prjcode_set, doc=r"""prjcode""")
    latbase = property(_wcs.WorldCoor_latbase_get, _wcs.WorldCoor_latbase_set, doc=r"""latbase""")
    ncoeff1 = property(_wcs.WorldCoor_ncoeff1_get, _wcs.WorldCoor_ncoeff1_set, doc=r"""ncoeff1""")
    ncoeff2 = property(_wcs.WorldCoor_ncoeff2_get, _wcs.WorldCoor_ncoeff2_set, doc=r"""ncoeff2""")
    zpnp = property(_wcs.WorldCoor_zpnp_get, _wcs.WorldCoor_zpnp_set, doc=r"""zpnp""")
    changesys = property(_wcs.WorldCoor_changesys_get, _wcs.WorldCoor_changesys_set, doc=r"""changesys""")
    printsys = property(_wcs.WorldCoor_printsys_get, _wcs.WorldCoor_printsys_set, doc=r"""printsys""")
    ndec = property(_wcs.WorldCoor_ndec_get, _wcs.WorldCoor_ndec_set, doc=r"""ndec""")
    degout = property(_wcs.WorldCoor_degout_get, _wcs.WorldCoor_degout_set, doc=r"""degout""")
    tabsys = property(_wcs.WorldCoor_tabsys_get, _wcs.WorldCoor_tabsys_set, doc=r"""tabsys""")
    rotmat = property(_wcs.WorldCoor_rotmat_get, _wcs.WorldCoor_rotmat_set, doc=r"""rotmat""")
    coorflip = property(_wcs.WorldCoor_coorflip_get, _wcs.WorldCoor_coorflip_set, doc=r"""coorflip""")
    offscl = property(_wcs.WorldCoor_offscl_get, _wcs.WorldCoor_offscl_set, doc=r"""offscl""")
    wcson = property(_wcs.WorldCoor_wcson_get, _wcs.WorldCoor_wcson_set, doc=r"""wcson""")
    naxis = property(_wcs.WorldCoor_naxis_get, _wcs.WorldCoor_naxis_set, doc=r"""naxis""")
    naxes = property(_wcs.WorldCoor_naxes_get, _wcs.WorldCoor_naxes_set, doc=r"""naxes""")
    wcsproj = property(_wcs.WorldCoor_wcsproj_get, _wcs.WorldCoor_wcsproj_set, doc=r"""wcsproj""")
    linmode = property(_wcs.WorldCoor_linmode_get, _wcs.WorldCoor_linmode_set, doc=r"""linmode""")
    detector = property(_wcs.WorldCoor_detector_get, _wcs.WorldCoor_detector_set, doc=r"""detector""")
    instrument = property(_wcs.WorldCoor_instrument_get, _wcs.WorldCoor_instrument_set, doc=r"""instrument""")
    ctype = property(_wcs.WorldCoor_ctype_get, _wcs.WorldCoor_ctype_set, doc=r"""ctype""")
    c1type = property(_wcs.WorldCoor_c1type_get, _wcs.WorldCoor_c1type_set, doc=r"""c1type""")
    c2type = property(_wcs.WorldCoor_c2type_get, _wcs.WorldCoor_c2type_set, doc=r"""c2type""")
    ptype = property(_wcs.WorldCoor_ptype_get, _wcs.WorldCoor_ptype_set, doc=r"""ptype""")
    units = property(_wcs.WorldCoor_units_get, _wcs.WorldCoor_units_set, doc=r"""units""")
    radecsys = property(_wcs.WorldCoor_radecsys_get, _wcs.WorldCoor_radecsys_set, doc=r"""radecsys""")
    radecout = property(_wcs.WorldCoor_radecout_get, _wcs.WorldCoor_radecout_set, doc=r"""radecout""")
    radecin = property(_wcs.WorldCoor_radecin_get, _wcs.WorldCoor_radecin_set, doc=r"""radecin""")
    eqin = property(_wcs.WorldCoor_eqin_get, _wcs.WorldCoor_eqin_set, doc=r"""eqin""")
    eqout = property(_wcs.WorldCoor_eqout_get, _wcs.WorldCoor_eqout_set, doc=r"""eqout""")
    sysin = property(_wcs.WorldCoor_sysin_get, _wcs.WorldCoor_sysin_set, doc=r"""sysin""")
    syswcs = property(_wcs.WorldCoor_syswcs_get, _wcs.WorldCoor_syswcs_set, doc=r"""syswcs""")
    sysout = property(_wcs.WorldCoor_sysout_get, _wcs.WorldCoor_sysout_set, doc=r"""sysout""")
    center = property(_wcs.WorldCoor_center_get, _wcs.WorldCoor_center_set, doc=r"""center""")
    wcsl = property(_wcs.WorldCoor_wcsl_get, _wcs.WorldCoor_wcsl_set, doc=r"""wcsl""")
    lin = property(_wcs.WorldCoor_lin_get, _wcs.WorldCoor_lin_set, doc=r"""lin""")
    cel = property(_wcs.WorldCoor_cel_get, _wcs.WorldCoor_cel_set, doc=r"""cel""")
    prj = property(_wcs.WorldCoor_prj_get, _wcs.WorldCoor_prj_set, doc=r"""prj""")
    lngcor = property(_wcs.WorldCoor_lngcor_get, _wcs.WorldCoor_lngcor_set, doc=r"""lngcor""")
    latcor = property(_wcs.WorldCoor_latcor_get, _wcs.WorldCoor_latcor_set, doc=r"""latcor""")
    distcode = property(_wcs.WorldCoor_distcode_get, _wcs.WorldCoor_distcode_set, doc=r"""distcode""")
    distort = property(_wcs.WorldCoor_distort_get, _wcs.WorldCoor_distort_set, doc=r"""distort""")
    command_format = property(_wcs.WorldCoor_command_format_get, _wcs.WorldCoor_command_format_set, doc=r"""command_format""")
    ltm = property(_wcs.WorldCoor_ltm_get, _wcs.WorldCoor_ltm_set, doc=r"""ltm""")
    ltv = property(_wcs.WorldCoor_ltv_get, _wcs.WorldCoor_ltv_set, doc=r"""ltv""")
    idpix = property(_wcs.WorldCoor_idpix_get, _wcs.WorldCoor_idpix_set, doc=r"""idpix""")
    ndpix = property(_wcs.WorldCoor_ndpix_get, _wcs.WorldCoor_ndpix_set, doc=r"""ndpix""")
    wcs = property(_wcs.WorldCoor_wcs_get, _wcs.WorldCoor_wcs_set, doc=r"""wcs""")
    wcsdep = property(_wcs.WorldCoor_wcsdep_get, _wcs.WorldCoor_wcsdep_set, doc=r"""wcsdep""")
    wcsname = property(_wcs.WorldCoor_wcsname_get, _wcs.WorldCoor_wcsname_set, doc=r"""wcsname""")
    wcschar = property(_wcs.WorldCoor_wcschar_get, _wcs.WorldCoor_wcschar_set, doc=r"""wcschar""")
    logwcs = property(_wcs.WorldCoor_logwcs_get, _wcs.WorldCoor_logwcs_set, doc=r"""logwcs""")

    def __init__(self):
        r"""__init__(self) -> WorldCoor"""
        _wcs.WorldCoor_swiginit(self, _wcs.new_WorldCoor())
    __swig_destroy__ = _wcs.delete_WorldCoor

# Register WorldCoor in _wcs:
_wcs.WorldCoor_swigregister(WorldCoor)



