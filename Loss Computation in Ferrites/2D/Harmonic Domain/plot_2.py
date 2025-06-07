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
P_exi = data["p_exi"]

f_mfem = data_sim["fc"]
P_eddy_mfem = data_sim['P_eddy']
P_exi_mfem = data_sim["P_mag"]

plt.figure()
plt.plot(f,P_eddy/1000,label="Peddy FreeFem")
plt.plot(f_mfem, P_eddy_mfem/1000, label="Peddy mfem")
plt.plot(f,P_exi/1000,label="Pexi FreeFem")
plt.plot(f_mfem, P_exi_mfem/1000, label="Pexi mfem")
plt.loglog()
plt.semilogx()
plt.grid()
plt.legend()

plt.show()


