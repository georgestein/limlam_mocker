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

### Plot results
llm.plot_results(mapinst,k,Pk,Pk_sampleerr,params)

llm.write_time('Finished Line Intensity Mapper')
