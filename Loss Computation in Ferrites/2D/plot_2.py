import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
import os

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, "build", "datafolder", "donnes_riccardo", "simul_N30_10mT.txt")
data_sim_path = os.path.join(current_path, "build", "datafolder", "data_tau_fixed.csv")

data = pd.read_csv(data_path, sep=" ")
data_sim = pd.read_csv(data_sim_path, sep=";")

f = data["f"]
P_eddy = data["p_eddy"]

f_mfem = data_sim["fc"]
P_eddy_mfem = data_sim['Ploss']

plt.figure()
plt.plot(f,P_eddy/1000,label="Peddy FreeFem")
plt.plot(f_mfem, P_eddy_mfem/2000, label="Peddy mfem")
# plt.loglog()
plt.semilogx()
plt.grid()
plt.legend()

plt.show()


