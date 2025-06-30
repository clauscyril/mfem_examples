import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_flux_1.csv")


height = 7.59e-3
w = 5.3e-3

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]/1000
flux = data["flux"]/w/height
phi = data['phi_imposed']/w/height
fluxH = data["fluxH"]/w/height
phiH = data['phiH_imposed']/w/height
NI = data['NI']

# p_eddy = np.array(p_eddy[500:])

print(np.mean(p_eddy[200:]))

height = 7.59e-3
w = 5.3e-3

plt.figure()
plt.plot(t, flux, label='Flux de B mesuré ', color='blue')
plt.plot(t, phi, label='Flux de B "imposé"', color='red', linestyle='--')
plt.legend()
plt.grid()
plt.figure()
plt.plot(t, fluxH, label='Flux de H mesuré', color='blue')
plt.plot(t, phiH, label='flux de H "imposé', color='red', linestyle='--')
plt.legend()
plt.grid()
# plt.figure()
# plt.plot(t, fluxH - phiH, label='Difference plux imposé et mesuré', color='blue')
# # plt.plot(t, phiH, label='flux de H "imposé', color='red', linestyle='--')
# plt.legend()
# plt.grid()
plt.figure()
plt.plot(t, NI, label='NI')
plt.legend()
plt.figure()
plt.plot(t, p_eddy,label="P_eddy")
plt.legend()

fig, ax1 = plt.subplots()
# Courbe 1 sur l'axe gauche
ax1.plot(t, flux, 'b-', label='B (normalized)')
ax1.set_xlabel('t')
ax1.set_ylabel('B', color='black')
ax1.tick_params(axis='y', labelcolor='b')

# Deuxième axe Y à droite
ax2 = ax1.twinx()
ax2.plot(t, fluxH, 'r--', label='H (normalized))')
ax2.set_ylabel('H', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Titre et affichage
plt.title("Double échelle Y")
plt.show()

plt.show()