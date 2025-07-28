import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


current_path = os.path.dirname(os.path.abspath(__file__))
path_h = os.path.join(current_path, "Harmonic Domain", "build", "datafolder", "data_tau_fixed.csv")
data_h = pd.read_csv(path_h, sep=";")

f_h = data_h['fc']/1000
p_h = data_h['P_eddy']/1000

f_t = np.array([100, 500, 800, 1000, 1800, 2000])
P_t = np.array([13.72065, 2307.70793,13528.6679, 29511.492, 152590.2899, 202404.767])/1000


plt.plot(f_t, P_t, "o", label="Time Domain")
plt.plot(f_h, p_h, "--", label="Harmonic Domain")
plt.semilogx()
plt.xlabel("F (kHz)")
plt.ylabel("P_eddy (KW/m3)")
# plt.loglog()
plt.grid()
plt.legend()
plt.show()


