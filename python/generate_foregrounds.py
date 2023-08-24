import litebird_sim as lbs
import sys
import healpy as hp
import numpy as np
import os
import matplotlib.pyplot as plt
from lbs2fgbuster import get_fgbuster_instrument_from_imo
import logging
logging.disable(logging.INFO)

def gen_input_maps(toml_filename, base_path):
    """
    generate FG, fg models listed in `fg_models`
    produced components:
     - Temperature : synch, dust, freefree, ame
     - Polarization: synch, dust

    Maps are produced in galactic, in K_CMB units
    Maps are saved in basepath+'/d0s0' directory
    """
    inst = get_fgbuster_instrument_from_imo()
    sim = lbs.Simulation(parameter_file=toml_filename, random_seed=33)

    nside       = int(sim.parameters["general"]["nside"])
    fg_synch    =     sim.parameters['fg_models']['fg_synch']
    fg_dust     =     sim.parameters['fg_models']['fg_dust']
    fg_freefree =     sim.parameters['fg_models']['fg_freefree']
    fg_ame      =     sim.parameters['fg_models']['fg_ame']

    Mbsparams   = lbs.MbsParameters( # for Polarization
                make_cmb         = False,
                make_fg          = True,
                fg_models        = [fg_synch, fg_dust], #freefree & ame are included later, but only for T
                gaussian_smooth  = False,
                bandpass_int     = False,
                maps_in_ecliptic = False,
                parallel_mc      = False,
                nside            = nside,
                units            = "K_CMB",
            )
    freq_channels = inst["channel"] #all freq channels including bandwidth
    telescope     = inst["telescope"]
    channel_id    = np.arange(len(freq_channels))
    
    for (i, ifreq, teles) in zip(channel_id, freq_channels, telescope):
        print("Channel: ", ifreq)
        mbs = lbs.Mbs(simulation=sim, 
                        parameters=Mbsparams, 
                        channel_list  = [lbs.FreqChannelInfo.from_imo(
                            sim.imo,
                            "/releases/v1.3/satellite/{}/{}/channel_info".format(teles, ifreq))]
                       )
        all_maps = mbs.run_all()[0]
        
        for j in range(3):
            if type(all_maps[ifreq][j,0])!=np.float64:
                print("Map doesn't have np.float64")
                break

        if not os.path.exists(base_path):
            os.mkdir(base_path)
        if not os.path.exists(base_path+'/d0s0'):
            os.mkdir(base_path+'/d0s0')
        if not os.path.exists(base_path+'/d0s0' + '/nside_'+str(nside)):
            os.mkdir(base_path+'/d0s0'+'/nside_'+str(nside))
        
        components = '_'.join(Mbsparams.fg_models)
        hp.write_map(
            filename  = base_path+'/d0s0'+'/nside_'+str(nside)+'/'+ifreq+"_"+components+"_"+'nside_'+str(nside)+'.fits',
            m         = all_maps[ifreq],
            dtype     = np.float64,
            overwrite = True,
            )
        del all_maps


#----------------------------------- MAIN -----------------------------------#
nside       = int(sys.argv[1])
npix        = hp.nside2npix(nside)
inst        = get_fgbuster_instrument_from_imo()
print("nside : ", nside)

fg_synch    = 'pysm_synch_0'    #COMPLETE HERE
fg_dust     = 'pysm_dust_0'     #COMPLETE HERE
fg_freefree = 'pysm_freefree_1' #COMPLETE HERE
fg_ame      = 'pysm_ame_1'
#base_path   = '/group/cmb/litebird/usr/ytakase/maps/foregrounds' #COMPLETE HERE where maps are saved
base_path   = './'

toml_filename = os.getcwd() + '/' + 'generate_inputs_nside_{}.toml'.format(nside)
with open(toml_filename, 'w') as f:
    f.write('[general]\n') # ---------- GENERAL ----------
    f.write('nside = '+str(nside)+'\n')
    f.write('[fg_models]\n') # ---------- FG_MODELS ----------
    f.write('fg_synch = \''+fg_synch+'\'\n')
    f.write('fg_dust = \''+fg_dust+'\'\n')
    f.write('fg_freefree = \''+fg_freefree+'\'\n')
    f.write('fg_ame = \''+fg_ame+'\'\n')
    f.close()

gen_input_maps(toml_filename, base_path)
os.remove(toml_filename)