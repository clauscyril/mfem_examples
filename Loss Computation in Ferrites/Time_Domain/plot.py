import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_1.csv")

height = 7.59e-3
w = 5.3e-3
print(10e-3*w*height)

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]/1000
p_exc = data['p_exc']/1000
flux = data["flux"]
NI = data['NI']
fluxH = data['fluxH']
phiH = data["phiH"]

# f = 100e3
# T = 1/f
# Ts = T/1000

# t = t[500:]
# p_eddy = np.array(p_eddy[500:])


print(np.mean(p_eddy[400:]) + np.mean(p_exc[400:]))    


# plt.figure()
# plt.plot(freq, np.abs(fft_p), label="Fréquentiel")
# plt.grid()

# plt.figure()
# plt.plot(t,p_eddy, label="Puissance instantanée (kW/m^3)")
# plt.xlabel('Time (s)')
# plt.ylabel("Power (kW/m3)")
# plt.grid()
# plt.legend()
# plt.figure()
# plt.plot(t,flux, label="Magnetic flux (Wb)")
# plt.plot(t, NI)
# plt.xlabel('Time (s)')
# plt.ylabel('Magnetic flux (Wb)')
# plt.legend()
# plt.grid()
# plt.show()


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(16, 12))

# Premier graphe
ax1.plot(t, flux, label='Magnetic Flux', color='blue')
ax1.set_title("Magnetic Flux")
ax1.set_ylabel("Magnetic Flux (Wb)")
# ax1.set_xlabel("Time (s)")
ax1.grid(True)
# ax1.legend()

# Deuxième graphe
ax2.plot(t, NI, label='I', color='red', linestyle='--')
ax2.set_title("NI (400 kHz)")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("NI (A)")
ax2.grid(True)
# ax2.legend()

# Ajustement automatique
plt.tight_layout()

plt.figure()
plt.plot(t,p_eddy, label="Power eddy")
plt.plot(t,p_exc, label="Power exc")
plt.grid()
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Power (kW/m3)")
plt.figure()
plt.plot(t,fluxH,  label="Flux de H")
plt.plot(t,phiH, '--', label="Flux de H from phi de b")
plt.grid()
plt.legend()
plt.show()