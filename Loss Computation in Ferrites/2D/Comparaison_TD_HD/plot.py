import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))

path_td = os.path.join(path, "build", "data")
path_hd = os.path.join(path, "build", "data", "HD_data.csv")

data = pd.read_csv(path_hd, sep=';')
f = data['f']
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


flux_mod_val = flux_mod[0]
angle_val = angle[0]
f_val = f[0] 
flux_hd_func = flux_mod_val * np.sqrt(2) *  np.sin(2*np.pi*f_val*t + angle_val)


print(np.size(p_eddy))
print(np.mean(p_eddy[np.size(p_eddy)//2:]))

plt.figure()
plt.plot(f, p_eddy_mean, label="Time Domain")
plt.plot(f, p_eddy_hd, label="Steady-sate")
plt.grid()
plt.legend()
plt.xlabel('f (kHz)')
plt.ylabel('Power (kW/m3)')

plt.figure()
# plt.plot(t, flux,'-',label="Time Domain")
plt.plot(t, flux_hd_func,'--', label="Harmonic steady-state")
plt.ylabel("Magnetic Flux (Wb)")
plt.xlabel("Time (s)")
# plt.title("Comparison of the magnetic flux in Time Domain and Harmonic steady-state for a sin wave of " + str(f_val/1000)+'kHz')
plt.legend(loc='upper right')
plt.show()