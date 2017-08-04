#!/usr/bin/env python
from __future__ import division
import numpy              as np
import matplotlib.pylab   as plt
import scipy              as sp
import limlam_mocker      as llm
#Get Parameters for run
import params             as params

llm.debug.verbose = True
llm.write_time('Starting Line Intensity Mapper')

### Setup maps to output
mapinst   = llm.params_to_mapinst(params);

### Load halos from catalogue
halos, cosmo = llm.load_peakpatch_catalogue(params.halo_catalogue_file)
halos        = llm.cull_peakpatch_catalogue(halos, params.min_mass, mapinst)

### Calculate Luminosity of each halo
halos.Lco    = llm.Mhalo_to_Lco(halos, params.model, params.coeffs)

### Bin halo luminosities into map
mapinst.maps = llm.Lco_to_map(halos,mapinst)

### Output map to file
llm.save_maps(mapinst)

### Calculate power spectrum
k,Pk,Nmodes = llm.map_to_pspec(mapinst,cosmo)
Pk_sampleerr = Pk/np.sqrt(Nmodes)


llm.write_time('Finished Line Intensity Mapper')


### Plot central frequency map
plt.rcParams['font.size'] = 16
if params.plot_cube:
    plt.figure()
    im = plt.imshow(np.log10(mapinst.maps[:,:,params.nmaps//2]+1e-6), extent=[-mapinst.fov_x/2,mapinst.fov_x/2,-mapinst.fov_y/2,mapinst.fov_y/2],vmin=-1,vmax=2)
    plt.colorbar(im,label=r'$log_{10}\ T_b\ [\mu K]$')
    plt.xlabel('degrees',fontsize=20)
    plt.ylabel('degrees',fontsize=20)
    plt.title('simulated map at {0:.3f} GHz'.format(mapinst.nu_bincents[params.nmaps//2]),fontsize=24)

if params.plot_pspec:
    plt.figure()
    plt.errorbar(k,k**3*Pk/(2*np.pi**2),k**3*Pk_sampleerr/(2*np.pi**2),
                    lw=3,capsize=0)
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.grid(True)
    plt.xlabel('k [1/Mpc]',fontsize=18)
    plt.ylabel('$\\Delta^2(k)$ [$\\mu$K$^2$]',fontsize=18)
    plt.title('simulated line power spectrum',fontsize=24)

if params.plot_cube or params.plot_pspec:
    plt.show()


