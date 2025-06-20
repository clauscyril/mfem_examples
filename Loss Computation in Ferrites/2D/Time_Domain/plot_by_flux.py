import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_flux_1.csv")

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]/1000
flux = data["flux"]
phi = data['phi_imposed']
fluxH = data["fluxH"]
phiH = data['phiH_imposed']
NI = data['NI']

# p_eddy = np.array(p_eddy[500:])

print(np.mean(p_eddy[200:]))



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
plt.figure()
plt.plot(t, fluxH - phiH, label='Flux de H mesuré', color='blue')
# plt.plot(t, phiH, label='flux de H "imposé', color='red', linestyle='--')
plt.legend()
plt.grid()
plt.figure()
plt.plot(t, NI)
plt.show()