import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))

path_td = os.path.join(path, "build", "data")
path_hd = os.path.join(path, "build", "data", "HD_data.csv")

data = pd.read_csv(path_hd, sep=';')
f = np.sort(np.array(data['f']))
flux_hd = data['flux_r'] + 1j* data["flux_i"]
angle = np.angle(flux_hd)
flux_mod = np.abs(flux_hd)
p_eddy_hd = data['p_eddy']

p_eddy_mean = []


for filename in os.listdir(path_td):
    if "TD_" in filename:
        full_path = os.path.join(path_td, filename)
        print(full_path)

        data_td = pd.read_csv(full_path, sep=';')
        t = data_td['t']
        p_eddy = data_td['p_eddy']
        flux = data_td['flux']
        p_mean = np.mean(p_eddy[np.size(p_eddy)//2:])
        p_eddy_mean.append(p_mean)
        # print(np.mean(p_eddy[np.size(p_eddy)//2:]))
        # p_eddy_mean.append(np.mean(p_eddy[np.size(p_eddy)//2:]))
        # np.append(p_eddy_mean, np.mean(p_eddy[np.size(p_eddy)//2:]))   

p_eddy_mean.sort()


flux_mod_val = flux_mod[len(flux_mod)-1]
angle_val = angle[len(flux_mod)-1]
f_val = f[len(flux_mod)-1] 




print(np.size(p_eddy))
print(np.mean(p_eddy[np.size(p_eddy)//2:]))

plt.figure()
plt.plot(f, p_eddy_mean, label="Time Domain")
plt.plot(f, p_eddy_hd, label="Steady-sate")
plt.grid()
plt.legend()
plt.xlabel('f (kHz)')
plt.ylabel('Power (kW/m3)')


full_path = os.path.join(path_td, "TD_1999.csv")
data_td = pd.read_csv(full_path, sep=';')
t = data_td['t']
p_eddy = data_td['p_eddy']
flux = data_td['flux']
flux_hd_func = flux_mod_val * np.sqrt(2) *  np.sin(2*np.pi*f_val*t + angle_val)

plt.figure()
# plt.plot(t, flux,'-',label="Time Domain")
plt.plot(t, flux_hd_func,'--', label="Harmonic steady-state")
plt.plot(t, flux, '-', label="Time Domain")
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel("Magnetic Flux (Wb)", fontsize=16)
plt.xlabel("Time (s)", fontsize=16)
plt.title("Magnetic flux, f = 2MHz", fontsize=16)
plt.legend(loc='upper right', fontsize=14)
plt.show()