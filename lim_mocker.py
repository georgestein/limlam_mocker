#!/usr/bin/env python
from __future__ import division
from src.mods     import *


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


### Plot first frequency map
if params.plot_cube:
    im = plt.imshow(np.log10(map.maps[:,:,params.nmaps//2]+1e-1), extent=[-map.fov_x/2,map.fov_x/2,-map.fov_y/2,map.fov_y/2])
    plt.colorbar(im,label=r'$log_{10}\ T_b\ [\mu K]$')
    plt.xlabel('degrees',fontsize=20)
    plt.ylabel('degrees',fontsize=20)
    plt.show() 
