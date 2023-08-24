def deconvolution(maps, fwhm, cut_off=191):
    """
    Deconvolve a set of input maps using a Gaussian beam.

    This function takes input maps, a Full Width at Half Maximum (FWHM) value
    for the Gaussian beam, and an optional cut-off multipole value. It performs
    deconvolution on the input maps using the specified Gaussian beam, up to the
    given cut-off multipole if provided.

    Parameters:
    -----------
        maps (array-like)      : Input maps to be deconvolved. If maps.shape is (3, npix),
                                 it's assumed to be a set of polarized maps for each Stokes parameter.
                                 Otherwise, it's assumed to be an temperature map.
        fwhm (float)           : Full Width at Half Maximum (FWHM) value of the Gaussian beam in radians.
        cut_off (int, optional): Cut-off multipole value. Multipole values above this cut-off
                                 will not be deconvolved. Default is 191 which was used in PTEP in the likelihood function.

    Returns:
    --------
        array-like             : Deconvolved maps after applying the Gaussian beam correction.
                                 The shape of the output maps will match the shape of the input maps.

    Note:
    -----
        The function internally uses the HEALPix library for spherical harmonic transformations.
    """
    nside = hp.get_nside(maps)
    npix  = hp.nside2npix(nside)
    lmax  = 3*nside - 1
    alm   = hp.map2alm(maps, use_weights=True)
    if maps.shape == (3, npix):
        bl         = hp.gauss_beam(fwhm=fwhm,lmax=lmax,pol=True)
        alm_deconv = np.zeros(alm.shape, dtype=complex)
        for m in range(lmax+1):
            for i in range(lmax-m+1):
                l = i + m
                if l <= cut_off:
                    idx = hp.Alm.getidx(lmax, l, m)
                    for stokes in range(3):
                        alm_deconv[stokes, idx] = alm[stokes,idx]/bl[l,stokes]
                else: 
                    for stokes in range(3):
                        alm_deconv[stokes, idx] = 0.0
    else:
        bl    = hp.gauss_beam(fwhm=fwhm,lmax=lmax,pol=False)
        alm_deconv = np.zeros(len(alm), dtype=complex)
        for m in range(lmax+1):
            idx1 = hp.Alm.getidx(lmax, m, m)
            idx2 = hp.Alm.getidx(lmax, lmax, m)
            alm_deconv[idx1:idx2+1] = alm[idx1:idx2+1] / bl[m:lmax+1]
        
    return hp.alm2map(alm_deconv, nside)