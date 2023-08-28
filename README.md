# lbsim_tools
## Installation 
- Dependency: [litebird_sim](https://github.com/litebird/litebird_sim)

If you haven't installed the litebird_sim you must install it before installing lbsim_tools.
```
$ git clone https://github.com/yusuke-takase/lbsim_tools.git
$ cd lbsim_tools
$ (lbs_env)$ python install -e .
```
## API
- `deconvolution(maps, fwhm, cut_off)`
    - This function provide the deconvolution for the map which is convolved by a Gaussian beam. The [usage](./notebooks/deconv_verification.ipynb) is availble with the verification of the function. 
- `get_fgbuster_instrument_from_imo(imo_version)`
    - This function genarates a table which is used for FGBuster by using the litebird_sim imo. 
- `c2d(cl, ell_start=2.)`
    - Convert $C_\ell$ to $D_\ell$.
- `get_planck_cmap()`
    - Generate planck color scheme.
        ```
        # Usage
        import lbsim_tools as lt
        import healpy as hp
        import numpy as np
        m    = np.arange(hp.nside2npix(32))
        cmap = lt.get_planck_cmap()
        hp.mollview(m, cmap=cmap)
        ```
- `read_fiducial_cl(r)`
    - This function reads the power spectrum of the CMB used in the map base simulation of litebird_sim. 
    It refers to the power spectrum calculated with the specified tensor-to-scalar ratio, $r$ by setting the argument `r` to `r=0` or `r=1`.
- `forecast(lmax, cl_sys, rmin=1e-8, rmax=1e-1, rresol=1e5, iter=0, verbose=False, test=False, bias=1e-5)`
    - This function estimates the tensor scalar ratio from the power spectrum using the likelihood function used in [PTEP: P88, Sec. (5.3.2)](https://academic.oup.com/ptep/article/2023/4/042F01/6835420). 
    In doing so, it excludes multipoles above the $\ell$ specified by the argument `lmax`. Enter the $B$-mode power spectrum of the systematic error in the argument `cl_sys` (The unit of it must be $\mu K_{CMB}^2$). 
    For example, the power spectrum of the map obtained from the difference between the input map and the output map including systematic effects i.e. residual map corresponds to this.
    



## Scripts 
The python files which are included in [script](./script) is executable from your terminal.
- detsfile_generator.py
    ```
    (lbs_env)$ python detsfile_generator.py
    ```
    - The FPU will be displaied that you requested by using the IMo.
    - You can check the information of the detector that you specified by clicking.
    - You can save the text file that contains the detectors list what you chose in e2e simulation format. 

- generate_foregrounds.py
    - Save foreground maps by using litebird_sim 
 
 