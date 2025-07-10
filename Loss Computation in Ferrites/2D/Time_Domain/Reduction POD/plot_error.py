import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

steps = np.linspace(0,490,50)

path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data")


errors = []

for el in steps :
    step = str(int(el))
    path_fom = os.path.join(path, "fom", "Hn_"+step+".csv")
    path_rom = os.path.join(path, "reduced", "Hn_reduced_"+step+"0")

    data_fom = pd.read_csv(path_fom)
    Hn_fom = data_fom['Hn']

    data_rom = pd.read_csv(path_rom, header=None)
    Hn_rom = data_rom[0]
    diff = Hn_fom - Hn_rom
    error = np.linalg.norm(diff)
    error_ref = np.linalg.norm(Hn_fom)
    errors.append(error/error_ref)


path_fom = os.path.join(path, "fom", "TD_1.csv")
path_rom = os.path.join(path, "reduced", "TD_1.csv")

data_fom = pd.read_csv(path_fom, sep=";")
data_rom = pd.read_csv(path_rom, sep=";")

print(data_rom)

flux_fom = data_fom["flux"]
flux_rom = data_rom["flux"]
flux_rom_reduced = data_rom['flux_reduced']
diff_flux = flux_fom - flux_rom

p_eddy_fom = data_fom["p_eddy"]
p_eddy_rom = data_rom["p_eddy"]

t = data_rom['t']


plt.plot(errors)
plt.title("Relative Error between Fom and Rom")
plt.xlabel("Iteration")
plt.ylabel("Relative error")
plt.grid()

plt.figure()

plt.plot(t, flux_rom, label="Flux Rom")
plt.plot(t, flux_rom_reduced, label="Flux Rom reduced")
# plt.plot(t, flux_fom,"--", label="Flux FOM")
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

plt.plot(t, p_eddy_rom, label="Power Rom")
# plt.plot(t, p_eddy_fom, label="Power FOM")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Eddy Current losses (W/m3))")
plt.legend()



plt.show()
