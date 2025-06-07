import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "build", "power_2000.csv")

data = pd.read_csv(path, sep=";")
t = data["t"]
p_eddy = data["p_eddy"]
flux = data["flux"]

# f = 100e3
# T = 1/f
# Ts = T/1000

# t = t[500:]
# p_eddy = np.array(p_eddy[500:])
print(len(t))
print(np.mean(p_eddy))


# plt.figure()
# plt.plot(freq, np.abs(fft_p), label="Fréquentiel")
# plt.grid()

plt.figure()
plt.plot(t,p_eddy, label="Simulation")
# plt.plot(t, envelope)
plt.figure()
plt.plot(t,flux)
plt.show()
