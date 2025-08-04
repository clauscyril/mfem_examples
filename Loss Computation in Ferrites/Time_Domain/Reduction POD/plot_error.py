import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data")

path_fom = os.path.join(path, "fom", "TD_0.csv")
path_rom = os.path.join(path, "reduced", "TD_1.csv")

data_fom = pd.read_csv(path_fom, sep=";")
data_rom = pd.read_csv(path_rom, sep=";")

flux_fom = data_fom["flux"]
flux_rom = data_rom["flux"]
diff_flux = flux_fom - flux_rom

p_eddy_fom = data_fom["p_eddy"]
p_eddy_rom = data_rom["p_eddy"] 

power_error = abs(p_eddy_rom - p_eddy_fom)

height = 7.59e-3
w = 5.3e-3

t = data_rom['t']
NI = data_fom["NI"]

print(f"Power Fom :  {np.mean(p_eddy_fom[200:])}")
print(f"Power ROM reduced : {np.mean(p_eddy_rom[200:])}")

plt.figure()

# plt.plot(t, flux_rom, label="Flux Rom")
plt.plot(t, flux_rom/w/height, label="Flux Rom reduced")
plt.plot(t, flux_fom/w/height,"--", label="Flux FOM")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Magnetic flux (Wb)")
plt.legend()

plt.figure()
plt.plot(t, diff_flux, label="Difference Flux (FOM - ROM)")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Magnetic flux (Wb)")
plt.legend()


plt.figure()
plt.plot(t, NI, label="NI")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("NI (A)")
plt.legend()


plt.figure()
# plt.plot(t, p_eddy_rom, 'black', label="Power Rom")
plt.plot(t, p_eddy_rom, 'green', label="Power Rom Reduced ")
plt.plot(t, p_eddy_fom, 'r-.', label="Power FOM")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Eddy Current losses (W/m3))")
plt.legend()


plt.figure()
plt.plot(t, power_error , label="Power difference")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Eddy Current losses (W/m3))")
plt.legend()



plt.show()
