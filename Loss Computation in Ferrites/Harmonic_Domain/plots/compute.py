import os
import subprocess
import pandas as pd
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path_prog = os.path.join(path,"build")
exec_path = os.path.join(path_prog, "simHD")

mesh_path = "../../mesh/square.msh"

f_list = np.linspace(100e3, 2e6, 19)
P_eddy = []
P_mag = []
P_tot = []
flux_real = []
flux_imag = []
NI = []

for f in f_list: 
    exec = subprocess.run([exec_path, "-p", mesh_path, "-f", str(f), "-mat", "N30", "-b", "10"], capture_output=True, cwd=path_prog, check=True)

    data = str(exec.stdout.strip()).replace("'", "").split('\\n')[-1].split(";")

    P_eddy.append(float(data[1]))
    P_mag.append(float(data[2]))
    P_tot.append(float(data[3]))
    flux_real.append(float(data[4]))
    flux_imag.append(float(data[5]))
    NI.append(float(data[6]))

plt.figure()
plt.plot(f_list, P_tot)
plt.grid()
plt.xlabel('f(Hz)')
plt.ylabel('Total losses (W/m3)')
plt.show()
