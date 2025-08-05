import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

path = os.path.dirname(os.path.abspath(__file__))
path_file_sim = os.path.join(path,"data", "MFEM", "Princeton_mfemN30_15.csv")
path_file_meas = os.path.join(path, "data","Princeton","N30","N30-Sinusoidal_Phi_15.csv")

data_sim = pd.read_csv(path_file_sim, sep=";")
data_meas = pd.read_csv(path_file_meas, sep=",")

f = data_sim["fc"]
print(data_meas)
Power_meas = data_meas['Power_Loss']
Power_sim = data_sim['P_tot']

plt.plot(f, Power_meas, label="Measurements")
plt.plot(f, Power_sim, '-*', label='MFEM')

plt.legend()
plt.grid()
plt.show()