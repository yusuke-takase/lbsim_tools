import litebird_sim as lbs
import numpy as np
import healpy as hp
import astropy
import matplotlib.pyplot as plt
import os
import astropy.units as u

def get_det_xy(det_info):
    q = det_info.quat
    r = np.array([0.,0.,1.,0.])
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]])
    if telescope == "LFT":
        q_z = lbs.quat_rotation_z(np.deg2rad(180))
        q_z_conj = np.array([-q_z[0], -q_z[1], -q_z[2], q_z[3]])
        lbs.quat_left_multiply(r, *q)
        lbs.quat_right_multiply(r, *q_conj)
        lbs.quat_left_multiply(r, *q_z)
        lbs.quat_right_multiply(r, *q_z_conj)
    else:
        lbs.quat_left_multiply(r, *q)
        lbs.quat_right_multiply(r, *q_conj)
    theta, phi = hp.vec2ang(r[0:3])
    x = np.rad2deg(theta) * np.cos(phi)
    y = np.rad2deg(theta) * np.sin(phi)
    return (x,y)

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

def store_dets_info(detector):
    global selected_detector_list
    selected_detector_list.append(detector)

def remove_dets_info(detector):
    global selected_detector_list
    if detector in selected_detector_list:
        selected_detector_list.remove(detector)

def on_plot_click(event):
    global selected_detector_list
    white = "#1f77b4"
    red   = "#b41f44"
    if not event.inaxes:
        return
    for i, (x, y) in enumerate(zip(X, Y)):
        distance = ((x - event.xdata)**2 + (y - event.ydata)**2)**0.5
        if distance < 0.5:
            detector = total_detector_list[i]
            if scatter[i].get_markerfacecolor() == white:
                scatter[i].set_markerfacecolor(red)
                scatter[i].set_markeredgecolor(red)
                scatter[i].set_marker("*")
                scatter[i].set_markersize(12)
                info_text = gen_info_text(detector)
                info_box.set_text(info_text)
                selected_detector_list.append(detector)
            elif scatter[i].get_markerfacecolor() == red:
                scatter[i].set_markerfacecolor(white)
                scatter[i].set_markeredgecolor(white)
                scatter[i].set_marker("o")
                scatter[i].set_markersize(8)
                info_text = gen_info_text(detector)
                info_box.set_text(info_text)
                if detector in selected_detector_list:
                    selected_detector_list.remove(detector)
            fig.canvas.draw()

def ask_yes_or_no():
    print("Do you want to make a detector list file? [y/n]")
    while True:
        ans = input(">>> ").lower()  # 入力を小文字に変換して統一
        if ans == 'y':
            # 'y' が入力された場合の処理
            print("Creating a detector list file...")
            # ここにファイル作成の処理を追加
            break  # ループを抜ける
        elif ans == 'n':
            # 'n' が入力された場合の処理
            print("No detector list file will be created.")
            break  # ループを抜ける
        else:
            # 'y' または 'n' 以外の入力があった場合の処理
            print("Invalid input. Please enter 'y' or 'n'.")
            print("Do you want to make a detector list file? [y/n]")
    return ans


imo         = lbs.Imo()
sim         = lbs.Simulation(random_seed=None)
imo_version = "v1.3"

print("Input telescope name: ")
telescope     = input(">>> ")
channel_list  = []
inst_info     = sim.imo.query("/releases/"+imo_version+"/satellite/"+telescope+"/instrument_info")
channel_list  = inst_info.metadata["channel_names"]

print("The availavle channels are: ", channel_list)
print("Input the channel name: ")
channel       = input(">>> ")

print("Do you make a detector list file? [y/n]")
ans = ask_yes_or_no()
if ans == "y":
    print("Specify the directory to save: ")
    base_path = input(">>> ")


channel_info  = lbs.FreqChannelInfo.from_imo(
	imo=imo,
	url='/releases/'+imo_version+'/satellite/'+telescope+'/'+channel+'/channel_info',
)

total_detector_list = []
for detector_name in channel_info.detector_names:
	total_detector_list.append(lbs.DetectorInfo.from_imo(
			imo=imo,
			url='/releases/'+imo_version+'/satellite/'+telescope+'/'+channel+'/'+detector_name+'/detector_info',
		)
	)
total_number_of_dets = len(total_detector_list)

selected_detector_list = []
X   = np.zeros(len(total_detector_list))
Y   = np.zeros(len(total_detector_list))

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(1, 2, 1, aspect="equal")
scatter = []
for i,detector in enumerate(total_detector_list):
    X[i], Y[i] = get_det_xy(detector)
    scatter.append(ax1.plot(X[i], Y[i], "o", markersize=8, color="#1f77b4")[0])
ax1.set_xlabel(r"$\theta\cos(\phi)$ [degrees]")
ax1.set_ylabel(r"$\theta\sin(\phi)$ [degrees]")
if telescope == "LFT":
    ax1.set_xlim(-10,10)
    ax1.set_ylim(-6,6)
else:
    ax1.set_xlim(-15,15)
    ax1.set_ylim(-15,15)

ax2 = fig.add_subplot(1, 2, 2, aspect="equal")
info_box = ax2.text(0.02, 0.98, "", transform=ax2.transAxes, va="top", ha="left", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
ax2.set_axis_off()

fig.canvas.mpl_connect("button_press_event", on_plot_click)

plt.tight_layout()
plt.show()

if ans == "y":
    filename = "detectors_"+telescope+"_"+channel+"_T+B.txt"
    selected_detector_list = sorted(selected_detector_list, key=lambda detector: detector.name)
    generate_text_file(filename, selected_detector_list)
