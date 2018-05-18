Limlam Mocker: A package for Line Intensity Mocks
=================================================

Limlam mocker is a barebones code to create line intensity maps from a given halo catalogue. 
It only requires `NumPy <http://www.numpy.org/>`_ and `SciPy <http://www.scipy.org/>`_ - no special astronomy packages.

To Run
------

all parameters are in ./params.py and are very self explanatory. Set them up as you wish and run with 

        >>> ./lim_mocker.py


This will load in the halo catalogue, assign luminosities to each halo, and bin them up into a 3D intensity map data cube of size (npix_x,npix_y,nmaps). This map will then be saved using the npz format - eg https://docs.scipy.org/doc/numpy-1.12.0/reference/generated/numpy.savez.html . This file contains all required info for the maps - field of view, pixel size, frequency of maps, etc...

A small sample halo catalogue is included, but many more of larger sizes are available by contacting the author. 

the basic workflow
------------------
The example script provided in `lim_mocker.py` shows the basic workflow:

- `import params` initialises the parameters for the mock map.
- `load_peakpatch_catalogue` loads the halos from the lightcone specified.
- `cull_peakpatch_catalogue` implements a mass cutoff.
- `Mhalo_to_Lco` calculates luminosities for all halos.
- `params_to_mapinst` generates the map instance from given parameters.
- `Lco_to_map` then populates this map with temperatures.
- `save_maps` saves the maps to the npz file as described above.
- `map_to_pspec` then also calculates a 3D spherically averaged power spectrum for this map.

To add your own L_CO(M,z,...) function
--------------------------------------
add it to halos_to_luminosity.py, following the ones already there eg:    

        >>> dict = {'Li':          Mhalo_to_Lco_Li,
        >>>        'Padmanabhan': Mhalo_to_Lco_Padmanabhan}
            
        >>> def Mhalo_to_Lco_Li(halos, scatter):
        >>>        ...

To add your own halo catalogue
------------------------------
use the template in load_halos.py

Advanced usage: using `limlam_mocker` directly in your own code
---------------------------------------------------------------
This repo is a complete program in itself, designed so that you can painlessly run `lim_mocker.py` and not have to worry about the guts. If you want to use code other than `lim_mocker.py` with this library, the `limlam_mocker` folder inside this repository will function as a self-contained Python package.

Depending on the OS and the Python version, what may work is to simply place this folder in your local site-packages directory or create a symlink there. So if the contents of e.g. `~/.local/lib/python3.5/site-packages/limlam_mocker` (the local site-packages directory for Python 3.5 in Linux) are equal to the contents of the `limlam_mocker` folder in this repo (not the root contents of this repo!), you should be able to take advantage of `limlam_mocker` functions from any Python script on your computer, regardless of where it is placed.

This code was written by George Stein    - gstein@cita.utoronto.ca
    with many additions by Dongwoo Chung - dongwooc@stanford.edu

A version with many more options and functions (useful power spectrum calculation too!) can be found at https://github.com/dongwooc/imapper2, written by Tony Li and Dongwoo Chung.

License
-------

Baseband is licensed under the GNU General Public License v3.0 - see the
``LICENSE`` file.

