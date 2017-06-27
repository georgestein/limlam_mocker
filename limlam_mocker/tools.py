from __future__ import print_function
import datetime
import numpy as np

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

def params_to_mapinst(params):
    map             = empty_table() # creates empty class to put map info into 

    map.output_file = params.map_output_file 

    map.nmaps  = int(params.nmaps)
    map.fov_y  = float(params.fov_x)
    map.fov_x  = float(params.fov_y)
    map.npix_x = int(params.npix_x)
    map.npix_y = int(params.npix_y)
    map.nu_i   = float(params.nu_i)
    map.nu_f   = float(params.nu_f)

    # get arrays describing the final intensity map to be output
    # map sky angle dimension
    map.pix_size_x      = map.fov_x/map.npix_x 
    map.pix_size_y      = map.fov_y/map.npix_y

    map.Ompix = (map.pix_size_x*np.pi/180)*(map.pix_size_y*np.pi/180) # pixel size to convert to brightness temp 

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
