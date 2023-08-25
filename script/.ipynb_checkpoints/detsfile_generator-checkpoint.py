import litebird_sim as lbs
import numpy as np
import healpy as hp
import astropy
import matplotlib.pyplot as plt
import os
import astropy.units as u

def get_det_xy(det_info):
    q = det_info.quat
    r = np.array([0., 0., 1., 0.])
    if telescope == "LFT":
        # Rotate the FPU 180 degrees for projection if LFT is selected
        q_z = lbs.quat_rotation_z(np.deg2rad(180))
        for q_val in [q, q_z]:
            q_conj = np.array([-q_val[0], -q_val[1], -q_val[2], q_val[3]])
            lbs.quat_left_multiply(r, *q_val)
            lbs.quat_right_multiply(r, *q_conj)
    else:
        q_conj = np.array([-q[0], -q[1], -q[2], q[3]])
        lbs.quat_left_multiply(r, *q)
        lbs.quat_right_multiply(r, *q_conj)
    theta, phi = hp.vec2ang(r[:3])
    x = np.rad2deg(theta) * np.cos(phi)
    y = np.rad2deg(theta) * np.sin(phi)
    return x, y

def gen_info_text(detector):
    info_text = """
    Detector info.
        name: {detector.name}
        wafer: {detector.wafer}
        pixel: {detector.pixel}
        pixtype: {detector.pixtype}
        channel: {detector.channel}
        sampling_rate_hz: {detector.sampling_rate_hz}
        fwhm_arcmin: {detector.fwhm_arcmin}
        ellipticity: {detector.ellipticity}
        bandcenter_ghz: {detector.bandcenter_ghz}
        net_ukrts: {detector.net_ukrts}
        pol_sensitivity_ukarcmin: {detector.pol_sensitivity_ukarcmin}
        fknee_mhz: {detector.fknee_mhz}
        fmin_hz: {detector.fmin_hz}
        alpha: {detector.alpha}
        pol: {detector.pol}
        orient: {detector.orient}
        quat: {detector.quat}
    """
    info_text = info_text.format(**locals())
    return info_text

def generate_text_file(filename, selected_detector_list):
    header = "# Telescope     Channel         IMO_NET         Number_det      Scaled_NET      Detector_name\n"
    selected_number_of_dets = len(selected_detector_list)
    assumed_mission_duration_yr = 1.
    scaling_factor = np.sqrt(assumed_mission_duration_yr/3. * selected_number_of_dets/total_number_of_dets)
    with open(base_path+"/"+filename, "w") as file:
        file.write(header)
        for detector in selected_detector_list:
            scaled_net = np.round(detector.net_ukrts * scaling_factor, 2)
            line = f"{telescope}\t\t{channel}\t\t{detector.net_ukrts}\t\t{selected_number_of_dets}/{total_number_of_dets}\t\t{scaled_net}\t\t{detector.name}\n"
            file.write(line)
    print('The "'+filename+'" is generated.')
    print("Location: ", base_path+"/"+filename)

def on_plot_click(event):
    if not event.inaxes:
        return
    else:
        global selected_detector_list
        global scatter
        blue  = "#1f77b4"
        red   = "#b41f44"
        distance = ((X - event.xdata)**2 + (Y - event.ydata)**2)**0.5
        if np.min(distance) < 0.3:
            sorted_indices = np.argsort(distance)
            indices = [sorted_indices[0], sorted_indices[1]]
            for idx in indices:
                detector = total_detector_list[idx]
                if scatter[idx].get_markerfacecolor() == blue:
                    scatter[idx].set_markerfacecolor(red)
                    scatter[idx].set_markeredgecolor(red)
                    scatter[idx].set_marker("*")
                    scatter[idx].set_markersize(12)
                    info_text = gen_info_text(detector)
                    info_box.set_text(info_text)
                    selected_detector_list.append(detector)
                elif scatter[idx].get_markerfacecolor() == red:
                    scatter[idx].set_markerfacecolor(blue)
                    scatter[idx].set_markeredgecolor(blue)
                    scatter[idx].set_marker("o")
                    scatter[idx].set_markersize(8)
                    info_text = gen_info_text(detector)
                    info_box.set_text(info_text)
                    if detector in selected_detector_list:
                        selected_detector_list.remove(detector)
                fig.canvas.draw()
                #print("NumOfdets: ", len(selected_detector_list))

def ask_yes_or_no():
    print("Do you want to make a detector list file? [y/n]")
    while True:
        ans = input(">>> ").lower()
        if ans == 'y':
            print("Create a detector list file.")
            break
        elif ans == 'n':
            print("No detector list file will be created.")
            break
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
            print("Do you want to make a detector list file? [y/n]")
    return ans


#---------------------- Main function ----------------------------#
imo         = lbs.Imo()
sim         = lbs.Simulation(random_seed=None)
imo_version = "v1.3"

print("Input telescope name (LFT, MFT or HFT): ")
telescope     = input(">>> ")
inst_info     = sim.imo.query("/releases/"+imo_version+"/satellite/"+telescope+"/instrument_info")
channel_list  = inst_info.metadata["channel_names"]

print("The availavle channels are: ", channel_list)
print("Input the channel name: ")
channel       = input(">>> ")

print("Do you make a detector list file? [y/n]")
ans           = ask_yes_or_no()
if ans == "y":
    print("Specify the directory to save: ")
    base_path = input(">>> ")
    if base_path.endswith('/'):
        base_path = base_path[:-1]

channel_info  = lbs.FreqChannelInfo.from_imo(
	imo = imo,
	url = '/releases/'+imo_version+'/satellite/'+telescope+'/'+channel+'/channel_info',
)

total_detector_list = []
for detector_name in channel_info.detector_names:
	total_detector_list.append(lbs.DetectorInfo.from_imo(
			imo = imo,
			url = '/releases/'+imo_version+'/satellite/'+telescope+'/'+channel+'/'+detector_name+'/detector_info',
		)
	)
total_number_of_dets = len(total_detector_list)

selected_detector_list = []
X       = np.zeros(len(total_detector_list))
Y       = np.zeros(len(total_detector_list))


# Make a figure
fig     = plt.figure(figsize=(10,8))
ax1     = fig.add_subplot(1, 2, 1, aspect="equal")
scatter = []
for i,detector in enumerate(total_detector_list):
    X[i], Y[i] = get_det_xy(detector)
    scatter.append(ax1.plot(X[i], Y[i], "o", markersize=8, color="#1f77b4")[0])
ax1.set_xlabel(r"$\theta\cos(\phi)$ [degrees]")
ax1.set_ylabel(r"$\theta\sin(\phi)$ [degrees]")
"""
if telescope == "LFT":
    ax1.set_xlim(-10,10)
    ax1.set_ylim(-6,6)
else:
    ax1.set_xlim(-15,15)
    ax1.set_ylim(-15,15)
"""
ax2 = fig.add_subplot(1, 2, 2, aspect="equal")
info_box = ax2.text(0.02, 0.98, "", transform=ax2.transAxes, va="top", ha="left", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
ax2.set_axis_off()

print("Click is available...")
fig.canvas.mpl_connect("button_press_event", on_plot_click)
#fig.canvas.mpl_connect("button_press_event", get_detector_index)
plt.tight_layout()
plt.show()

# Save the detector list file.
if ans == "y":
    filename = "detectors_"+telescope+"_"+channel+"_T+B.txt"
    selected_detector_list = sorted(selected_detector_list, key=lambda detector: detector.name)
    generate_text_file(filename, selected_detector_list)
