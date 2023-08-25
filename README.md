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
 
 