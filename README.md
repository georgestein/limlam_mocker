This is as simple as possible code to create line intensity maps from a given halo catalogue. it only requires python, numpy, and scipy - no special astronomy packages.

TO RUN:
all parameters are in ./params.py and are very self explanatory. Set them up as you wish and run with 
$ ./lim_mocker.py

This will load in the halo catalogue, assign luminosities to each halo, and bin them up into a 3D Intensity map data cube of size (npix_x,npix_y,nmaps). This map will then be saved using the npz format - eg https://docs.scipy.org/doc/numpy-1.12.0/reference/generated/numpy.savez.html . This file contains all required into for the maps - field of view, pixel size, frequency of maps, etc...

A small sample halo catalogue is included, but many more of larger sizes are available by contacting the authour. 

TO ADD YOUR OWN LCO(M,z,...) FUNCTION:
add it to src/halos_to_luminosity.py, following the example:    
if model=='Li':
        return Mhalo_to_Lco_Li(halos)

TO ADD YOUR OWN HALO CATALOGUE:
copy src/load_halos.py

This code was written by George Stein    - gstein@cita.utoronto.ca
And with many additions by Dongwoo Chung - dongwooc@stanford.edu

A version with many more options and functions (useful powerspectrum calculatiosn too!) can be found at https://github.com/dongwooc/imapper2, written by Tony Li and Dongwoo Chung.

