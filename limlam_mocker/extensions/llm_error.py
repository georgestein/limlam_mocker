import numpy as np

def pspec_err_helper(mapinst,Tsys,Nfeeds,tobs,fwhm,cosmo,Wbvec=False):
    # Tsys in K; tobs in sec
    # Oobs in sr; fwhm, dpix in rad
    # dnu (channel BW), Dnu (total BW), nu_rest all in GHz
    z = mapinst.nu_rest/np.mean(mapinst.nu_binedges)-1
    dk = np.mean(np.diff(mapinst.k)) # 1/Mpc
    ctd = cosmo.comoving_transverse_distance(z).value
    dx_fwhm = ctd*fwhm/2.355 # Mpc
    dz = 299792.458/cosmo.H(z).value*mapinst.dnu/mapinst.nu_rest*(1+z)**2 # Mpc
    Oobs = mapinst.fov_x*mapinst.fov_y*(np.pi/180)**2
    signsq = 1e3*Tsys**2/(Nfeeds*mapinst.dnu*tobs*mapinst.Ompix/Oobs)
    Pn = signsq*ctd**2*mapinst.Ompix*dz
    kx,ky,kz = mapinst.kvec; kperpsq = kx**2+ky**2
    W_integrand = np.exp(-kperpsq*dx_fwhm**2)
    kbins, kgrid = mapinst.kbins, mapinst.kgrid
    Wbeam = np.histogram(kgrid[kgrid>0],bins=kbins,weights=W_integrand[kgrid>0])[0]/mapinst.Nmodes
    if Wbvec:
        return Pn,Wbeam,W_integrand
    else:
        return Pn,Wbeam

def snr_linespec(Pk,Pn,Nmodes,W):
    return np.sqrt(np.sum((Pk/(Pk+Pn)*np.sqrt(Nmodes)*W)**2))

def snr_xspec(Px,Pk,Pn,Pg,ngal_mean,Nmodes,Wx):
    return np.sqrt(np.sum(((Px*Wx)**2/(
                    ((Pk+Pn)*(Pg+1/ngal_mean)+Px**2)/2/Nmodes))))

def snr_rofk(Px,Pk,Pn,Pg,ngal_mean,Nmodes,W,Wx):
    rk = Px/(Pk*Pg)**0.5
    sigrk = rk/np.sqrt(Nmodes)*np.sqrt(
                (1+(Pk+Pn)*(Pg+1/ngal_mean)/Px**2)/2/Wx**2
                    + (1+Pn/Pk)**2/(4*W**2) + (1+1/ngal_mean/Pg)**2/4)
    return np.sqrt(np.sum((rk/sigrk)**2))
