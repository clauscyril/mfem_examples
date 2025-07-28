import os
import subprocess
import pandas as pd
import scipy.io
import numpy as np
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


exec_cpp_code = False

if exec_cpp_code:
    ### Executing c++ compiled program to compute the simulation
    path_prog = os.path.join(path,"build")
    exec_path = os.path.join(path_prog, "main")

    mats = ["N30", "T38"]
    b_peaks = ["10", "20"]

    # Parameters for the simulation 
    for mat in mats:
        for b_peak in b_peaks :
            subprocess.run([exec_path, "-m", mat, "-b", b_peak], cwd=path_prog, check=True) 

### Ploting the results

## Get measurments data
meas_N30 = scipy.io.loadmat(os.path.join(path,"data", "UNIPG", 'meas_N30_jan_2025.mat'))
meas_T38 = scipy.io.loadmat(os.path.join(path,"data", "UNIPG", 'meas_T38_jan_2025.mat'))
# meas_N87 = scipy.io.loadmat(os.path.join(path,"data", "UNIPG", 'meas_N87_jan_2025.mat'))

structN30 = meas_N30['meas']
dataN30 = structN30['data'][0, 0]
f = dataN30['f']
Losses = dataN30['loss']
f_N30_10 = np.array([e[0] for e in f[0,0]])
P_N30_10 = np.array([e[0] for e in Losses[0,0]])/1000
f_N30_20 = np.array([e[0] for e in f[0,1]])
P_N30_20 = np.array([e[0] for e in Losses[0,1]])/1000

structT38 = meas_T38['meas']
dataT38 = structT38['data'][0, 0]
f = dataT38['f']
Losses = dataT38['loss']
f_T38_10 = np.array([e[0] for e in f[0,0]])
P_T38_10 = np.array([e[0] for e in Losses[0,0]])/1000
f_T38_20 = np.array([e[0] for e in f[0,1]])
P_T38_20 = np.array([e[0] for e in Losses[0,1]])/1000

## Get simulation data
sim_N30_10 = pd.read_csv(os.path.join(path, "data", "MFEM", "unipg_mfemN30_10.csv"), sep=";")
sim_f_N30_10 = sim_N30_10['fc']
sim_P_N30_10 = sim_N30_10['P_tot']/1000

sim_N30_20 = pd.read_csv(os.path.join(path, "data", "MFEM", "unipg_mfemN30_20.csv"), sep=";")
sim_f_N30_20 = sim_N30_20['fc']
sim_P_N30_20 = sim_N30_20['P_tot']/1000

## Get simulation data
sim_T38_10 = pd.read_csv(os.path.join(path, "data", "MFEM", "unipg_mfemT38_10.csv"), sep=";")
sim_f_T38_10 = sim_T38_10['fc']
sim_P_T38_10 = sim_T38_10['P_tot']/1000

sim_T38_20 = pd.read_csv(os.path.join(path, "data", "MFEM", "unipg_mfemT38_20.csv"), sep=";")
sim_f_T38_20 = sim_T38_20['fc']
sim_P_T38_20 = sim_T38_20['P_tot']/1000



plt.figure()
plt.plot(f_N30_10, P_N30_10, label="Measurements")
plt.plot(sim_f_N30_10, sim_P_N30_10, label='MFEM')
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel('Power losses (kW/m3)')
plt.title("N30, b_peak = 10 mT")
plt.grid()

plt.figure()
plt.plot(f_N30_20, P_N30_20, label="Measurements")
plt.plot(sim_f_N30_20, sim_P_N30_20, label='MFEM')
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel('Power losses (kW/m3)')
plt.title("N30, b_peak = 20 mT")
plt.grid()


plt.figure()
plt.plot(f_T38_10, P_T38_10, label="Measurements")
plt.plot(sim_f_T38_10, sim_P_T38_10, label='MFEM')
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel('Power losses (kW/m3)')
plt.title("T38, b_peak = 10 mT")
plt.grid()

plt.figure()
plt.plot(f_T38_20, P_T38_20, label="Measurements")
plt.plot(sim_f_T38_20, sim_P_T38_20, label='MFEM')
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel('Power losses (kW/m3)')
plt.title("T38, b_peak = 20 mT")
plt.grid()



plt.show()