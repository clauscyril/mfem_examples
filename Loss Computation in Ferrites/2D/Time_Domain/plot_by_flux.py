import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_flux_1.csv")

data = pd.read_csv(path, sep=";")
t = data["t"][:15]
p_eddy = data["p_eddy"][:15]
flux = data["flux"][:15]
phi = data['phi_imposed'][:15]
fluxH = data["fluxH"][:15]
phiH = data['phiH_imposed'][:15]
NI = data['NI'][:15]

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
plt.figure()
plt.plot(t, p_eddy)
plt.show()