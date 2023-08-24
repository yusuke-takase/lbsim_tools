import pandas as pd
import litebird_sim as lbs

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
