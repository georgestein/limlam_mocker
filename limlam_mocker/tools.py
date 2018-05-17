from __future__ import print_function
from . import debug
import time
import datetime
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

class empty_table():
    """
    simple Class creating an empty table
    used for halo catalogue and map instances
    """
    def __init__(self):
        pass

    def copy(self):
        """@brief Creates a copy of the table."""
        return copy.copy(self)

def write_time(string_in):
    """
    write time info in as nicely formatted string
    """ 
    fmt       = '%H:%M:%S on %m/%d/%Y'
    timestamp = datetime.datetime.now().strftime(fmt)
    bar = 72*'-'
    print( '\n\n'+bar )
    print( string_in )
    print( 'Time:      '+timestamp )
    print( bar+'\n' )

    return

def timeme(method):
    """
    writes the time it takes to run a function
    To use, pput above a function definition. eg:
    @timeme
    def Lco_to_map(halos,map):
    """
    def wrapper(*args, **kw):
        startTime = int(round(time.time()))
        result = method(*args, **kw)
        endTime = int(round(time.time()))
                      
        print('  ',endTime - startTime,'sec')
        return result

    return wrapper

def params_to_mapinst(params):
    """
    Adds input parameters to be kept by the map class and gets map details

    Returns
    -------
    map : class
       contains all information about the map that the halos will be binned into
    """ 
    map             = empty_table() # creates empty class to put map info into 

    map.output_file = params.map_output_file 

    map.nmaps   = int(params.nmaps)
    map.fov_y   = float(params.fov_x)
    map.fov_x   = float(params.fov_y)
    map.npix_x  = int(params.npix_x)
    map.npix_y  = int(params.npix_y)
    map.nu_i    = float(params.nu_i)
    map.nu_f    = float(params.nu_f)
    map.nu_rest = float(params.nu_rest)
    map.z_i     = map.nu_rest/map.nu_i - 1
    map.z_f     = map.nu_rest/map.nu_f - 1

    # get arrays describing the final intensity map to be output
    # map sky angle dimension
    map.pix_size_x = map.fov_x/map.npix_x 
    map.pix_size_y = map.fov_y/map.npix_y

    # pixel size to convert to brightness temp 
    map.Ompix = (map.pix_size_x*np.pi/180)*(map.pix_size_y*np.pi/180) 

    map.pix_binedges_x = np.arange(-map.fov_x/2,map.fov_x/2+ map.pix_size_x, map.pix_size_x)
    map.pix_binedges_y = np.arange(-map.fov_y/2,map.fov_y/2+ map.pix_size_y, map.pix_size_y)

    map.pix_bincents_x =  0.5*(map.pix_binedges_x[:-1] + map.pix_binedges_x[:-1])
    map.pix_bincents_y =  0.5*(map.pix_binedges_y[:-1] + map.pix_binedges_y[:-1])

    # map frequency dimension 
    # negative steps as larger observed frequency means lower redshift
    map.dnu         = (map.nu_i - map.nu_f)/(map.nmaps)
    map.nu_binedges = np.arange(map.nu_i,map.nu_f-map.dnu,-map.dnu) 
    map.nu_bincents = map.nu_binedges[:-1] - map.dnu/2

    return map



# Cosmology Functions
# Explicitily defined here instead of using something like astropy 
# in order for ease of use on any machine 
def hubble(z,h,omegam):
    """
    H(z) in units of km/s
    """
    return h*100*np.sqrt(omegam*(1+z)**3+1-omegam)

def drdz(z,h,omegam):
    return 299792.458 / hubble(z,h,omegam)  

def chi_to_redshift(chi, cosmo):
    """
    Transform from redshift to comoving distance
    Agrees with NED cosmology to 0.01% - http://www.astro.ucla.edu/~wright/CosmoCalc.html
    """
    zinterp = np.linspace(0,4,10000)
    dz      = zinterp[1]-zinterp[0]

    chiinterp  = np.cumsum( drdz(zinterp,cosmo.h,cosmo.Omega_M) * dz)
    chiinterp -= chiinterp[0]
    z_of_chi   = sp.interpolate.interp1d(chiinterp,zinterp)

    return z_of_chi(chi)

def redshift_to_chi(z, cosmo):
    """
    Transform from comoving distance to redshift 
    Agrees with NED cosmology to 0.01% - http://www.astro.ucla.edu/~wright/CosmoCalc.html
    """
    zinterp = np.linspace(0,4,10000)
    dz      = zinterp[1]-zinterp[0]

    chiinterp  = np.cumsum( drdz(zinterp,cosmo.h,cosmo.Omega_M) * dz)
    chiinterp -= chiinterp[0]
    chi_of_z   = sp.interpolate.interp1d(zinterp,chiinterp)

    return chi_of_z(z)


def plot_results(mapinst,k,Pk,Pk_sampleerr,params):
    """
    Plot central frequency map and or powerspectrum
    """
    if debug.verbose: print("\n\tPlotting results")
    
    plt.rcParams['font.size'] = 14
    if params.plot_cube:
        plt.figure()
        im = plt.imshow(np.log10(mapinst.maps[:,:,params.nmaps//2]+1e-6), 
                        extent=[-mapinst.fov_x/2,mapinst.fov_x/2,-mapinst.fov_y/2,mapinst.fov_y/2],
                        vmin=-1,vmax=2)
        
        plt.colorbar(im,label=r'$log_{10}\ T_b\ [\mu K]$')
        plt.xlabel('degrees')
        plt.ylabel('degrees')
        plt.title('simulated map at {0:.3f} GHz'.format(mapinst.nu_bincents[params.nmaps//2]))

    if params.plot_pspec:
        plt.figure()
        plt.errorbar(k,k**3*Pk/(2*np.pi**2),k**3*Pk_sampleerr/(2*np.pi**2),
                     lw=3,capsize=0)
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')
        plt.grid(True)
        plt.xlabel('k [1/Mpc]')
        plt.ylabel('$\\Delta^2(k)$ [$\\mu$K$^2$]')
        plt.title('simulated line power spectrum')

    if params.plot_cube or params.plot_pspec:
        plt.show()

    return
