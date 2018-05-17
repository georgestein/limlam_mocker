from __future__ import absolute_import, print_function
import numpy as np
from . import debug
from .tools import *

@timeme
def map_to_pspec(map,cosmo):
    """
    Calculates power spectrum from 3D temperature cube

    Parameters
    ----------
    halos : class
        Contains all halo information (position, redshift, etc..)
    cosmo : class
        Contains all cosmology information (Omega_i, sigme_8, etc)   

    Returns
    -------
    k : `~numpy.ndarray`
        Wavevectors
    Pk : `~numpy.ndarray`
        Power at each k value
    nmodes : `~numpy.ndarray`
        Number of modes in each k bin
    """ 
    if debug.verbose: print("\n\tCalculating power spectrum")

    x,y,z = map.pix_binedges_x, map.pix_binedges_y, map.nu_binedges
    t     = map.maps

    zco   = redshift_to_chi(map.nu_rest/z-1,cosmo)
    # assume comoving transverse distance = comoving distance
    # (i.e. no curvature)

    avg_ctd = np.mean(zco) 
    xco     = x/(180)*np.pi*avg_ctd
    yco     = y/(180)*np.pi*avg_ctd

    dxco, dyco, dzco = [np.abs(np.mean(np.diff(d))) for d in (xco, yco, zco)]
    Pk_3D            = np.abs(np.fft.rfftn(t)*dxco*dyco*dzco)**2/np.abs(np.ptp(xco)*np.ptp(yco)*np.ptp(zco))

    kx        = 2*np.pi*np.fft.fftfreq(xco.size-1,d=dxco)
    ky        = 2*np.pi*np.fft.fftfreq(yco.size-1,d=dyco)
    kz        = 2*np.pi*np.fft.rfftfreq(zco.size-1,d=dzco)

    kgrid     = np.sqrt(sum(ki**2 for ki in np.meshgrid(kx,ky,kz,indexing='ij')))
    dk        = max(np.diff(kx)[0],np.diff(ky)[0],np.diff(kz)[0])

    kmax_dk   = int(np.ceil(max(np.amax(kx),np.amax(ky),np.amax(kz))/dk))
    kbins     = np.linspace(0,kmax_dk*dk,kmax_dk+1)

    Pk_nmodes = np.histogram(kgrid[kgrid>0],bins=kbins,weights=Pk_3D[kgrid>0])[0]
    nmodes    = np.histogram(kgrid[kgrid>0],bins=kbins)[0]

    Pk = Pk_nmodes/nmodes
    k  = (kbins[1:]+kbins[:-1])/2

    return k,Pk,nmodes
