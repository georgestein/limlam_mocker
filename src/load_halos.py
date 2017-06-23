import numpy as np
from   tools import *

def load_peakpatch_catalogue(filein):

    halos      = empty_table()            # creates empty class to put any halo info into  
    halo_info  = np.load(filein)     
    print "\thalo catalogue contains:\n\t\t", halo_info.files

    halos.M          = halo_info['M']     # halo mass in Msun
    
    halos.x_pos      = halo_info['x']     # halo position in comoving Mpc 
    halos.y_pos      = halo_info['y']
    halos.z_pos      = halo_info['z']
    halos.chi        = np.sqrt(halos.x_pos**2+halos.y_pos**2+halos.z_pos**2)

    halos.vx         = halo_info['vx']    # halo velocity in km/s
    halos.vy         = halo_info['vy']
    halos.vz         = halo_info['vz']

    halos.redshift   = halo_info['zhalo'] # redshift of halo
    halos.zformation = halo_info['zform'] # formation redshift of halo

    halos.nhalo = len(halos.M)
    
    halos.ra         = np.arctan(halos.x_pos/halos.z_pos)*180./np.pi
    halos.dec        = np.arctan(halos.y_pos/halos.z_pos)*180./np.pi

    print '\n\t%d halos loaded' % halos.nhalo

    return halos


def cull_peakpatch_catalogue(halos, min_mass):

    dm = [halos.M > min_mass]

    halos.M     = halos.M[dm]

    halos.x_pos = halos.x_pos[dm]
    halos.y_pos = halos.y_pos[dm]
    halos.z_pos = halos.z_pos[dm]
    halos.chi   = halos.chi[dm]

    halos.vx    = halos.vx[dm]
    halos.vy    = halos.vy[dm]
    halos.vz    = halos.vz[dm]

    halos.redshift   = halos.redshift[dm]
    halos.zformation = halos.zformation[dm]

    halos.ra  = halos.ra[dm]
    halos.dec = halos.dec[dm]

    halos.nhalo = len(halos.M)

    print '\n\t%d halos remain after mass cut' % halos.nhalo

    return halos
