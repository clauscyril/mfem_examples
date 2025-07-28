import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


path = os.path.dirname(os.path.abspath(__file__))

# Measurements from mat file
# Loss Computation in Ferrites\measurements\meas_N30_jan_2025.mat
mat_data = scipy.io.loadmat(os.path.join(path,"data", "UNIPG", 'meas_N30_jan_2025.mat'))
struct = mat_data['meas']
data = struct['data'][0, 0]
f = data['f']
Losses = data['loss']
f1 = np.array([e[0] for e in f[0,0]])/1000
Losses1 = np.array([e[0] for e in Losses[0,0]])/1000


path = os.path.dirname(os.path.abspath(__file__))

path_file = os.path.join(path, "data", "MFEM", "unipg_mfem.csv")

data = pd.read_csv(path_file, sep=";")

Fc = data['fc']/1000
Ploss = data['P_tot']/1000
P_eddy = data['P_eddy']/1000
P_mag = data['P_mag']/1000
flux = data['flux_r']

# Création du graphique
plt.figure()
plt.plot(f1, Losses1, label="Measurements")
plt.plot(Fc, P_eddy, '--.', label="eddy losses MFEM")
plt.plot(Fc, P_mag, '--.', label="exitation losses MFEM")
plt.plot(Fc, Ploss, '--.', label="Total losses MFEM")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')
# plt.yscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.title('Ferrite N87')
plt.grid()
plt.legend()


plt.figure()
plt.plot(Fc, flux, '--.', label="Tau 0")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Flux')
plt.title('Ferrite N87')
plt.grid()
plt.legend()
# Affichage du graphique

plt.show()