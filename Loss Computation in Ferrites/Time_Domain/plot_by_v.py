import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_v_1.csv")


height = 7.59e-3
w = 5.3e-3

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]/1000
flux = data["flux"]/w/height
V = data['V_imposed']
fluxH = data["fluxH"]/w/height
sigmaV = data['simgaV']
NI = data['NI']

# p_eddy = np.array(p_eddy[500:])

print(np.mean(p_eddy[int(len(p_eddy)/2):]))

height = 7.59e-3
w = 5.3e-3

plt.figure()
# plt.plot(t, flux, label='Flux de B mesuré ', color='blue')
plt.plot(t, V, label='V imposed', color='red', linestyle='--')
plt.legend()
plt.grid()
plt.figure()
# plt.plot(t, fluxH, label='Flux de H mesuré', color='blue')
plt.plot(t, sigmaV, label='sigma V', color='red', linestyle='--')
plt.legend()
plt.grid()
# plt.figure()
# plt.plot(t, fluxH - phiH, label='Difference plux imposé et mesuré', color='blue')
# # plt.plot(t, phiH, label='flux de H "imposé', color='red', linestyle='--')
# plt.legend()
# plt.grid()
plt.figure()
plt.plot(t, NI, label='NI')
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Current (A)")
plt.legend()
plt.figure()
plt.plot(t, p_eddy,label="P_eddy")
plt.grid()
plt.xlabel("Time (s)")
plt.ylabel("Power (kW/m3)")
plt.legend()

fig, ax1 = plt.subplots()
# Courbe 1 sur l'axe gauche
ax1.plot(t, flux, 'b-', label='B')
ax1.set_xlabel('t')
ax1.set_ylabel('B', color='black')
ax1.tick_params(axis='y', labelcolor='b')
plt.legend(loc='upper left')
# Deuxième axe Y à droite
ax2 = ax1.twinx()
ax2.plot(t, fluxH, 'r-', label='H')
ax2.set_ylabel('H', color='r')
ax2.tick_params(axis='y', labelcolor='r')
plt.legend()
# Titre et affichage
plt.title("H and B (normalized)")
plt.grid()
plt.show()

