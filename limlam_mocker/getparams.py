from __future__ import absolute_import
import params
from .pykDict import *

def getparams(filename):
    
    dict = pykDict()
    dict.read_from_file(filename)

    params.lco_model = dict['lco_model']
    params.scatter   = dict['scatter']
    
    params.halo_catalogue_file = dict['halo_catalogue_file']
    params.min_mass            = dict['min_mass']

    params.nu_rest = dict['nu_rest']
    params.nu_i    = dict['nu_i']
    params.nu_f    = dict['nu_f']

    params.nmaps  = dict['nmaps']
    params.fov_x  = dict['fov_x']
    params.fov_y  = dict['fov_y']
    params.npix_x = dict['npix_x']
    params.npix_y = dict['npix_y']

    params.map_output_file = dict['map_output_file']
    params.plot_cube       = dict['plot_cube']

    return
