from __future__ import absolute_import, print_function
import numpy as np
from  .tools import *
from . import debug

@timeme
def load_peakpatch_catalogue(filein):
    halos      = empty_table()            # creates empty class to put any halo info into  
    cosmo      = empty_table()            # creates empty class to put any cosmology info into  

    halo_info  = np.load(filein)     
    if debug.verbose: print("\thalo catalogue contains:\n\t\t", halo_info.files)
    
    #get cosmology from halo catalogue
    params_dict    = halo_info['cosmo_header'][()]
    cosmo.Omega_M  = params_dict.get('Omega_M')
    cosmo.Omega_B  = params_dict.get('Omega_B')
    cosmo.Omega_L  = params_dict.get('Omega_L')
    cosmo.h        = params_dict.get('h'      )
    cosmo.ns       = params_dict.get('ns'     )
    cosmo.sigma8   = params_dict.get('sigma8' )

    halos.M          = halo_info['M']     # halo mass in Msun
    
    halos.x_pos      = halo_info['x']     # halo position in comoving Mpc 
    halos.y_pos      = halo_info['y']
    halos.z_pos      = halo_info['z']
    halos.chi        = np.sqrt(halos.x_pos**2+halos.y_pos**2+halos.z_pos**2)

    halos.vx         = halo_info['vx']    # halo velocity in km/s
    halos.vy         = halo_info['vy']
    halos.vz         = halo_info['vz']

    halos.redshift   = halo_info['zhalo'] # observed redshift incl velocities
                                          
    halos.zformation = halo_info['zform'] # formation redshift of halo

    halos.nhalo = len(halos.M)
    
    halos.ra         = np.arctan2(-halos.y_pos,halos.z_pos)*180./np.pi
    halos.dec        = np.arcsin(  halos.x_pos/halos.chi  )*180./np.pi

    if debug.verbose: print('\n\t%d halos loaded' % halos.nhalo)

    return halos, cosmo

@timeme
def cull_peakpatch_catalogue(halos, min_mass, mapinst):
    dm = [(halos.M > min_mass) * (halos.redshift >= mapinst.z_i)
                                * (np.abs(halos.ra) <= mapinst.fov_x/2)
                                * (np.abs(halos.dec) <= mapinst.fov_y/2)
                                * (halos.redshift <= mapinst.z_f)]

    for i in dir(halos):
        if i[0]=='_': continue
        try:
            setattr(halos,i,getattr(halos,i)[dm])
        except TypeError:
            pass
    halos.nhalo = len(halos.M)

    if debug.verbose: print('\n\t%d halos remain after mass/map cut' % halos.nhalo)

    return halos
