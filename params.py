# This is the parameter file for all halo
# and map parameters

### Lco(M, z, ...) model
model   = 'Li' # models includes are "Li" and 'Padmanabhan'
coeffs  = None # specify None for default coeffs

### Halo parameters
halo_catalogue_file = 'catalogues/peakpatchcatalogue_1pt4deg_z2pt4-3pt4.npz'
min_mass            = 2.5e10

### Map parameters
nu_rest = 115.27 # rest frame frequency of CO(1-0) transition in GHz
nu_i    = 34.    # GHz
nu_f    = 26.

nmaps   = 100
fov_x   = 1.4    # in degrees
fov_y   = 1.4    # in degrees
npix_x  = 256
npix_y  = 256 

map_output_file = './Lco_cube_trial'

### Plot parameters 
plot_cube = True
plot_pspec = True
