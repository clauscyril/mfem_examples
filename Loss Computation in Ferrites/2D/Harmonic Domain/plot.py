import scipy.io
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os



path = os.path.dirname(os.path.abspath(__file__))

# Measurements from mat file
mat_data = scipy.io.loadmat(os.path.join(path, "..", "measurements", 'meas_N30_jan_2025.mat'))
struct = mat_data['meas']
data = struct['data'][0, 0]
f = data['f']
I = data['Imax']
Losses = data['loss']
f1 = np.array([e[0] for e in f[0,0]])/1000
I1 = np.array([e[0] for e in I[0,0]])
Losses1 = np.array([e[0] for e in Losses[0,0]])/1000








path = os.path.dirname(os.path.abspath(__file__))

path_file = os.path.join(path,"build", "datafolder", "data_tau_fixed.csv")
path_file0 = os.path.join(path,"build", "datafolder", "data0.csv")
path_file1 = os.path.join(path,"build", "datafolder", "data1.csv")
path_file2 = os.path.join(path,"build", "datafolder", "data2.csv")
path_file3 = os.path.join(path,"build", "datafolder", "data3.csv")
path_file4 = os.path.join(path,"build", "datafolder", "data4.csv")

data = pd.read_csv(path_file, sep=";")
data0 = pd.read_csv(path_file0, sep=";")
data1 = pd.read_csv(path_file1, sep=";")
data2 = pd.read_csv(path_file2, sep=";")
data3 = pd.read_csv(path_file3, sep=";")
data4 = pd.read_csv(path_file4, sep=";")



Fc = data['fc']/1000
Ploss = data['Ploss']/1000
flux = data['flux']
I = data['Imax']
I = I * I1[0]/I[0]


Fc0 = data0['fc']/1000
Ploss0 = data0['Ploss']/1000
flux0 = data0['flux']

Fc1 = data1['fc']/1000
Ploss1 = data1['Ploss']/1000
flux1 = data1['flux']

Fc2 = data2['fc']/1000
Ploss2 = data2['Ploss']/1000
flux2 = data2['flux']

Fc3 = data3['fc']/1000
Ploss3 = data3['Ploss']/1000
flux3 = data3['flux']

Fc4 = data4['fc']/1000
Ploss4 = data4['Ploss']/1000
flux4 = data4['flux']


# Création du graphique
plt.figure()
plt.plot(Fc, Ploss, '--.', label="Tau 0")
# plt.plot(f1, Losses1)
# plt.plot(Fc1, Ploss1, '--.', label="Tau 1")
# plt.plot(Fc2, Ploss2, '--.', label="Tau 2")
# plt.plot(Fc3, Ploss3, '--.', label="Tau 3")
# plt.plot(Fc4, Ploss4, '--.', label="Tau 4")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')
plt.yscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.title('Ferrite N87')
plt.grid()
plt.legend()


plt.figure()
plt.plot(Fc, flux, '--.', label="Tau 0")
# plt.plot(f1, I1)
# plt.plot(Fc1, flux1, '--.', label="Tau 1")
# plt.plot(Fc2, flux2, '--.', label="Tau 2")
# plt.plot(Fc3, flux3, '--.', label="Tau 3")
# plt.plot(Fc4, flux4, '--.', label="Tau 4")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Flux')
plt.title('Ferrite N87')
plt.grid()
plt.legend()
# Affichage du graphique



# plt.figure()
# plt.plot(Fc, flux, '--.', label="Tau 0")
# # plt.plot(Fc1, flux1, '--.', label="Tau 1")
# # plt.plot(Fc2, flux2, '--.', label="Tau 2")
# # plt.plot(Fc3, flux3, '--.', label="Tau 3")
# # plt.plot(Fc4, flux4, '--.', label="Tau 4")

# # Définition de l'échelle logarithmique pour l'axe y
# plt.xscale('log')

# # Ajout de labels et titre
# plt.xlabel('Frequency (kHz)')
# plt.ylabel('Flux')
# plt.title('Ferrite N87')
# plt.grid()
# plt.legend()
# # Affichage du graphique


plt.show()