import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.interpolate
import sys

def Mhalo_to_Lco(halos, model, scatter):

    dict = {'Li':          Mhalo_to_Lco_Li,
            'Padmanabhan': Mhalo_to_Lco_Padmanabhan}

    if model in dict.keys():
        return dict[model](halos, scatter)

    else:
        sys.exit('\n\n\tYour model, '+model+', does not seem to exist\n\t\tPlease check src/halos_to_luminosity.py to add it\n\n')


def Mhalo_to_Lco_Li(halos, scatter):
    """
    halo mass to SFR to L_CO 
    following the Tony li 2016 model
    arXiv 1503.08833
    """
    # Power law parameters from paper
    delta_mf  = 1.0
    alpha     = 1.37
    beta      = -1.74
    sigma_sfr = 0.3
    sigma_lco = 0.3

    # Get Star formation rate
    sfr_interp_tab = get_sfr_table()
    sfr            = sfr_interp_tab.ev(np.log10(halos.M), np.log10(halos.redshift+1))
    if scatter:
        sfr            = add_log_normal_scatter(sfr, sigma_sfr)

    # infrared luminosity
    lir      = sfr * 1e10 / delta_mf
    alphainv = 1./alpha
    # Lco' (observers units)
    Lcop     = lir**alphainv * 10**(-beta * alphainv)
    # Lco in L_sun
    Lco      =  4.9e-5 * Lcop
    if scatter:
        Lco      = add_log_normal_scatter(Lco, sigma_lco)

    print '\n\tMhalo to Lco calculated'

    return Lco 

def Mhalo_to_Lco_Padmanabhan(halos, scatter):
    """
    halo mass to L_CO 
    following the Padmanabhan 2017 model
    arXiv 1706.01471
    """
    m10,m11,n10,n11,b10,b11,y10,y11 = 4.17e12,-1.17,0.0033,0.04,0.95,0.48,0.66,-0.33
    z = halos.redshift
    hm = halos.M
    m1 = 10**(np.log10(m10)+m11*z/(z+1));
    n = n10+n11*z/(z+1);
    b = b10+b11*z/(z+1);
    y = y10+y11*z/(z+1);
    Lprime = 2*n*hm/((hm/m1)**(-b)+(hm/m1)**y);
    L = Lprime#*4.9e-5*(rest_freq/CO_REST_FREQ)**3;
    return L.flatten()

def get_sfr_table():
    """
    LOAD SFR TABLE 
    Columns are: z+1, logmass, logsfr, logstellarmass
    Intermediate processing of tabulated data      
    """

    dat_zp1, dat_logm, dat_logsfr, _ = np.loadtxt('tables/sfr_behroozi_release.dat',
                                                  unpack=True)

    dat_logzp1 = np.log10(dat_zp1)
    dat_sfr    = 10.**dat_logsfr

    # Reshape arrays
    dat_logzp1  = np.unique(dat_logzp1)    # log(z), 1D 
    dat_logm    = np.unique(dat_logm)    # log(Mhalo), 1D        
    dat_sfr     = np.reshape(dat_sfr, (dat_logm.size, dat_logzp1.size))

    # Get interpolated SFR value(s)    
    sfr_interp_tab = sp.interpolate.RectBivariateSpline(dat_logm, dat_logzp1, dat_sfr,
                                                        kx=1, ky=1)
    return sfr_interp_tab


def add_log_normal_scatter(data,dex):
    """
    Return array x, randomly scattered by a log-normal distribution with sigma=dexscatter. [via @tonyyli - https://github.com/dongwooc/imapper2]
    Note: scatter maintains mean in linear space (not log space).
    """
    if (dex<=0):
        return data
    # Calculate random scalings
    sigma       = dex * 2.302585 # Stdev in log space (DIFFERENT from stdev in linear space), note: ln(10)=2.302585
    mu          = -0.5*sigma**2

    randscaling = np.random.lognormal(mu, sigma, data.shape)
    xscattered  = np.where(data > 0, data*randscaling, data)

    return xscattered
