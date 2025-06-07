import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd


path = os.path.dirname(os.path.abspath(__file__))

# Measurements from mat file
mat_data = scipy.io.loadmat(os.path.join(path,'../../measurements', 'meas_N30_jan_2025.mat'))
struct = mat_data['meas']
data = struct['data'][0, 0]
f = data['f']
I = data['Imax']
Losses = data['loss']
f1 = np.array([e[0] for e in f[0,0]])/1000
I1 = np.array([e[0] for e in I[0,0]])
Losses1 = np.array([e[0] for e in Losses[0,0]])/1000

# Simulation results
print(path)
path_file = os.path.join(path,"build","datafolder", "data_tau_fixed.csv")
print(path_file)
data0 = pd.read_csv(path_file, sep=";")

F = data0['fc']/1000
P_eddy= data0['P_eddy']/1000
P_mag = data0["P_mag"]/1000
P_tot = data0["P_tot"]/1000
flux = data0['flux']

plt.figure()
plt.plot(f1, Losses1,'-o', label="Measurements")
plt.plot(F, P_tot,'-o', label="MFEM total losses 2D")
plt.plot(F, P_eddy,'--', label="MFEM Eddy")
plt.plot(F, P_mag,'--', label="MFEM mag")
plt.legend()
plt.xscale('log')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.grid()

plt.figure()
plt.plot(f1, Losses1,'-o', label="Measurements")
plt.plot(F, P_tot,'-o', label="MFEM total losses 2D")
plt.plot(F, P_eddy,'--', label="MFEM Eddy")
plt.plot(F, P_mag,'--', label="MFEM mag")
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.grid()
plt.show()




# np.savetxt(os.path.join(path,"build", 'f.csv'), f1, delimiter=';')
# print(I1)
# # Sauvegarder le tableau I dans un fichier CSV
# np.savetxt(os.path.join(path,"build", 'I.csv'), I1, delimiter=';')

