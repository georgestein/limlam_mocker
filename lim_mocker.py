#!/usr/bin/env python
from __future__ import division
from limlam_mocker.mods import *


### Load halos from catalogue
halos     = load_peakpatch_catalogue(params.halo_catalogue_file)
halos     = cull_peakpatch_catalogue(halos, params.min_mass)


### Calculate line freq from redshift
halos.nu  = params.nu_rest/(halos.redshift+1)       


### Calculate Luminosity of each halo
halos.Lco = Mhalo_to_Lco(halos, params.model, params.scatter)


### Bin halo luminosities into map
map.maps  = Lco_to_map(halos,map)

### Output map to file
save_maps(map)

### Calculate power spectrum
## if the catalogue contains halo comoving distances and redshifts ...
##     we don't need astropy!
from scipy.interpolate import UnivariateSpline
redshift_to_chi = UnivariateSpline(halos.redshift,halos.chi)
k,Pk,Pk_sampleerr = map_to_pspec(map,redshift_to_chi,params.nu_rest)

### Plot first frequency map
if params.plot_cube:
    plt.figure()
    im = plt.imshow(np.log10(map.maps[:,:,params.nmaps//2]+1e-1), extent=[-map.fov_x/2,map.fov_x/2,-map.fov_y/2,map.fov_y/2])
    plt.colorbar(im,label=r'$log_{10}\ T_b\ [\mu K]$')
    plt.xlabel('degrees',fontsize=20)
    plt.ylabel('degrees',fontsize=20)

if params.plot_pspec:
    plt.figure()
    plt.errorbar(k,k**3*Pk/(2*np.pi**2),k**3*Pk_sampleerr/(2*np.pi**2),
                    capsize=0)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.xlabel('k (1/Mpc)',fontsize=20)
    plt.ylabel('$\\Delta^2(k)$ ($\\mu$K$^2$)',fontsize=20)

if params.plot_cube or params.plot_pspec:
    plt.show()
