from __future__ import print_function
from . import debug
import time
import datetime
import numpy as np
import scipy as sp

class empty_table():
    """ 
    brief Class describing a table.    
    """
    def __init__(self):
        pass

    def copy(self):
        """ 
        @brief Creates a copy of the table.          
        """
        return copy.copy(self)

def write_time(string_in):
    fmt       = '%H:%M:%S on %m/%d/%Y'
    timestamp = datetime.datetime.now().strftime(fmt)
    bar = 72*'-'
    print( '\n\n'+bar )
    print( string_in )
    print( 'Time:      '+timestamp )
    print( bar+'\n' )

    return

def timeme(method):
    def wrapper(*args, **kw):
        startTime = int(round(time.time()))
        result = method(*args, **kw)
        endTime = int(round(time.time()))
                      
        if debug.verbose: print('  ',endTime - startTime,'sec')
        return result

    return wrapper

def params_to_mapinst(params):
    map             = empty_table() # creates empty class to put map info into 

    map.output_file = params.map_output_file 

    map.nmaps  = int(params.nmaps)
    map.fov_x  = float(params.fov_x)
    map.fov_y  = float(params.fov_y)
    map.npix_x = int(params.npix_x)
    map.npix_y = int(params.npix_y)
    map.nu_i   = float(params.nu_i)
    map.nu_f   = float(params.nu_f)
    map.nu_rest= float(params.nu_rest)
    map.z_i    = map.nu_rest/map.nu_i - 1
    map.z_f    = map.nu_rest/map.nu_f - 1

    # get arrays describing the final intensity map to be output
    # map sky angle dimension
    map.pix_size_x      = map.fov_x/map.npix_x 
    map.pix_size_y      = map.fov_y/map.npix_y

    map.Ompix = (map.pix_size_x*np.pi/180)*(map.pix_size_y*np.pi/180) # pixel size to convert to brightness temp 

    map.pix_binedges_x = np.linspace(-map.fov_x/2,map.fov_x/2,map.npix_x+1)
    map.pix_binedges_y = np.linspace(-map.fov_y/2,map.fov_y/2,map.npix_y+1)

    map.pix_bincents_x =  0.5*(map.pix_binedges_x[1:] + map.pix_binedges_x[:-1])
    map.pix_bincents_y =  0.5*(map.pix_binedges_y[1:] + map.pix_binedges_y[:-1])

    # map frequency dimension 
    # use linspace to ensure nmaps channels
    map.nu_binedges = np.linspace(map.nu_i,map.nu_f,map.nmaps+1) 
    map.dnu         = np.mean(np.diff(map.nu_binedges))
    map.nu_bincents = map.nu_binedges[:-1] - map.dnu/2
    return map



# Cosmology Functions
# Explicitily defined here instead of using something like astropy 
# in order for ease of use on any machine 

def hubble(z,h,omegam):
    return h*100*np.sqrt(omegam*(1+z)**3+1-omegam)

def drdz(z,h,omegam):
    return 299792.458 / hubble(z,h,omegam)  

def chi_to_redshift(chi, cosmo):
    # Transform from redshift to comoving distance
    # Agrees with NED cosmology to 0.01% - http://www.astro.ucla.edu/~wright/CosmoCalc.html
    zinterp = np.linspace(0,4,10000)
    dz      = zinterp[1]-zinterp[0]

    chiinterp  = np.cumsum( drdz(zinterp,cosmo.h,cosmo.Omega_M) * dz)
    chiinterp -= chiinterp[0]
    z_of_chi   = sp.interpolate.interp1d(chiinterp,zinterp)

    return z_of_chi(chi)

def redshift_to_chi(z, cosmo):
    # Transform from comoving distance to redshift 
    # Agrees with NED cosmology to 0.01% - http://www.astro.ucla.edu/~wright/CosmoCalc.html
    zinterp = np.linspace(0,4,10000)
    dz      = zinterp[1]-zinterp[0]

    chiinterp  = np.cumsum( drdz(zinterp,cosmo.h,cosmo.Omega_M) * dz)
    chiinterp -= chiinterp[0]
    chi_of_z   = sp.interpolate.interp1d(zinterp,chiinterp)

    return chi_of_z(z)


