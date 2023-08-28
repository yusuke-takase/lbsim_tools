import os 
from pathlib import Path
import numpy as np
import pandas as pd
import litebird_sim as lbs
import healpy as hp
from matplotlib.colors import ListedColormap

def get_fgbuster_instrument_from_imo(imo_version="v1.3"):
    """
    This function generates DataFrame which is used for FGBuster as `instrument` from IMo.
    """
    sim = lbs.Simulation(random_seed=None)
    telescopes     = ["LFT", "MFT", "HFT"]
    channel_list   = []
    freq           = [] 
    depth_p        = [] 
    fwhm           = []
    telescope_list = []
    
    for i in telescopes:
        inst_info = sim.imo.query("/releases/"+imo_version+"/satellite/"+i+"/instrument_info")
        channel_list.append(inst_info.metadata["channel_names"])
    channel_list = [item for sublist in channel_list for item in sublist]
    
    for i in channel_list:
        if i[0]   == "L":
            telescope = "LFT"
        elif i[0] == "M":
            telescope = "MFT"
        elif i[0] == "H":
            telescope = "HFT"
        chinfo = lbs.FreqChannelInfo.from_imo(sim.imo,
                  "/releases/{}/satellite/{}/{}/channel_info".format(imo_version, telescope, i))
        freq.append(chinfo.band.bandcenter_ghz)
        depth_p.append(chinfo.pol_sensitivity_channel_ukarcmin)
        fwhm.append(chinfo.fwhm_arcmin)
        telescope_list.append(telescope)
        
    instrument = pd.DataFrame(data = {
        'channel'    : channel_list,
        'frequency'  : freq,
        'depth_p'    : depth_p,
        'fwhm'       : fwhm,
        'f_sky'      : [1.0                  for i in range(len(channel_list))],
        'status'     : ["forecast"           for i in range(len(channel_list))],
        'reference'  : ["IMo-" + imo_version for i in range(len(channel_list))],
        'type'       : ["satellite"          for i in range(len(channel_list))],
        'experiment' : ["LiteBIRD"           for i in range(len(channel_list))],
        'telescope'  : telescope_list
    })
    #instrument = instrument.sort_values('frequency')
    return instrument

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


def get_planck_cmap():
    datautils_dir = Path(__file__).parent / "datautils"
    color_data = np.loadtxt(datautils_dir / "Planck_Parchment_RGB.txt")
    colombi1_cmap = ListedColormap(color_data/255.)
    colombi1_cmap.set_bad("gray")
    colombi1_cmap.set_under("white")
    planck_cmap = colombi1_cmap
    return planck_cmap

def c2d(cl, ell_start=2.):
    """ The function to convert C_ell to D_ell
    
    Parameters
    ----------
    cl: 1d-array
        Power spectrum
    ell_start:float (default = 2.)
        The multi-pole ell value of first index of the `cl`.
        
    Return
    ------
    dl: 1d-array
    """
    ell = np.arange(ell_start, len(cl)+ell_start)
    return cl*ell*(ell+1.)/(2.*np.pi)

def d2c(dl, ell_start=2.):
    """ The function to convert D_ell to C_ell
    
    Parameters
    ----------
    dl: 1d-array
        (Reduced) Power spectrum
    ell_start:float (default = 2.)
        The multi-pole ell value of first index of the `dl`.
        
    Return
    ------
    cl: 1d-array
    """
    ell = np.arange(ell_start, len(dl)+ell_start)
    return dl*(2.*np.pi)/(ell*(ell+1.))

def read_fiducial_cl(r=0):
    """ This function reads the power spectrum of the CMB used in the map base simulation of litebird_sim. 
    
    Parameter
    ---------
    r: int
    
    Return
    ------
    cl: 2d-arrays
    """
    datautils_dir = Path(lbs.__file__).parent / "datautils"
    if int(r) == 0:
        cl            = hp.read_cl(datautils_dir / "Cls_Planck2018_for_PTEP_2020_r0.fits")
    if int(r) == 1:
        cl            = hp.read_cl(datautils_dir / "Cls_Planck2018_for_PTEP_2020_tensor_r1.fits")
    return cl


def forecast(lmax, cl_sys, rmin=1e-8, rmax=1e-1, rresol=1e5, iter=0, verbose=False, test=False, bias=1e-5):
    """ The function to estimate the bias of tensor-to-scalar ratio.
    This function based on the PTEP paper: https://academic.oup.com/ptep/article/2023/4/042F01/6835420
    P88, Sec. (5.3.2)
    
    Usage and detail of the function
    --------------------------------
    The argument `rmin` and `rmax` represent a range for a first r survery. `rresol` is the resulution of the grid of r within the range.
    If the argument `iter` does not equal 0, the estimation is continued around the r which is estimated on a previous survey. 
    You can see the survey log with `verbose` makes True.
    If the argument `test=True` we can verify the correctness of this function. The estimation result should be same with the value we set as `bias`.

    
    Parameters
    ----------
    lmax   : int
    cl_sys : 1d-array
    rmin   : float
    rmax   : float
    rresol : float
    iter   : int
    verbose: bool
    test   : bool
    bias   : float
    
    Return
    ------
    data: Dict
    """
    rresol = int(rresol)
    gridOfr = np.linspace(rmin, rmax, num=rresol)
    # Load the dat file which includes model power spectrum
    #datautils_dir = Path(lbs.__file__).parent / "datautils"
    cl_r0 = read_fiducial_cl(r=0)
    cl_r1 = read_fiducial_cl(r=1)
    
    # Note that the dat file has a power spectrum value from ell = 2 to ~4000.
    # In order to keep the formalism, we insert zeros at ell=1,2 of power spectrum.
    cl_lens = cl_r0[2,:]
    cl_tens = cl_r1[2,:] 
    
    if test == True:
        print("The test option is True...")
        # cl_sys is replaced to cl_tens which is maltiplyed `bias`
        cl_sys[0:lmax+1] = cl_tens[0:lmax+1] * bias
        
    ell = np.arange(2, lmax+1)
    Nell = len(ell)
    delta_r = 0.
    likelihood = 0.
    grid_of_r_for_likelihood = 0.

    for j in range(iter + 1):
        Nr = len(gridOfr)
        likelihood = np.zeros(Nr)
        
        for i, grid_val in enumerate(gridOfr):
            Cl_hat = cl_sys[ell] + cl_lens[ell]
            Cl = grid_val * cl_tens[ell] + cl_lens[ell]
            likelihood[i] = np.sum((-0.5) * (2.*ell + 1.) * ((Cl_hat / Cl) + np.log(Cl) - ((2.*ell - 1.) / (2.*ell + 1.)) * np.log(Cl_hat)))
        
        likelihood = np.exp(likelihood - np.max(likelihood))
        maxid = np.argmax(likelihood)
        delta_r = gridOfr[maxid]
        survey_range = [delta_r - delta_r*(0.5/(j+1.)), delta_r + delta_r*(0.5/(j+1.))]
        gridOfr_old = gridOfr
        gridOfr = np.linspace(survey_range[0], survey_range[1], num=int(1e4))
        
        if verbose == True:
            print("*--------------------------- iter =", j, "---------------------------*")
            print("Δr                :", delta_r)
            print("Next survey range :", survey_range)
    
    # Calcurate the likelihood function again in the range that is delta_r*1e-3 < delta_r < delta_r*3.
    # Note that the delta_r has already been estimated, this likelihood is used for display. 
    grid_of_r_for_likelihood = np.linspace(delta_r*1e-3, delta_r*3., num=int(1e4))
    Nr = len(grid_of_r_for_likelihood)
    likelihood = np.zeros(Nr)
    
    for i, grid_val in enumerate(grid_of_r_for_likelihood):
        Cl_hat = cl_sys[ell] + cl_lens[ell]
        Cl = grid_val * cl_tens[ell] + cl_lens[ell]
        likelihood[i] = np.sum((-0.5) * (2.*ell + 1.) * ((Cl_hat / Cl) + np.log(Cl) - ((2.*ell - 1.) / (2.*ell + 1.)) * np.log(Cl_hat)))
    
    likelihood = np.exp(likelihood - np.max(likelihood))
    data = {"delta_r":delta_r, "grid_r":grid_of_r_for_likelihood, "likelihood":likelihood}
    return data
