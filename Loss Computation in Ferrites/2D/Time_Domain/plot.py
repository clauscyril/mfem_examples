import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "data", "TD_100.csv")

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]/1000
flux = data["flux"]

# f = 100e3
# T = 1/f
# Ts = T/1000

# t = t[500:]
# p_eddy = np.array(p_eddy[500:])
print(np.mean(p_eddy[500:]))


# plt.figure()
# plt.plot(freq, np.abs(fft_p), label="Fréquentiel")
# plt.grid()

plt.figure()
plt.plot(t,p_eddy, label="Puissance instantanée (kW/m^3)")
plt.xlabel('Time (s)')
plt.ylabel("Power (kW/m3)")
plt.grid()
plt.legend()
plt.figure()
plt.plot(t,flux, label="Magnetic flux (Wb)")
plt.xlabel('Time (s)')
plt.ylabel('Magnetic flux (Wb)')
plt.legend()
plt.grid()
plt.show()
