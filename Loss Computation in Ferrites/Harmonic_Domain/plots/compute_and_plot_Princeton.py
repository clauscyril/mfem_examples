import os
import subprocess
import pandas as pd
import scipy.io
import numpy as np
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_path = os.path.join(path, "data", "Princeton")

exec_cpp_code = True

if exec_cpp_code:
    ### Executing c++ compiled program to compute the simulation
    path_prog = os.path.join(path,"build")
    exec_path = os.path.join(path_prog, "sim")

    for dossier, sous_dossiers, fichiers in os.walk(data_path):
        for fichier in fichiers:
            mat = os.path.basename(fichier).split("-")[0]
            if mat == 'N30' :
                full_path = os.path.join(dossier, fichier)
                subprocess.run([exec_path, "-m", mat, "-p", full_path], cwd=path_prog, check=True) 

            

### Ploting the results

## Get measurments data

for dossier, sous_dossiers, fichiers in os.walk(data_path):
    for fichier in fichiers:
        mat = fichier.split('-')[0]
        phi = fichier.split('_')[-1].split('.')[0]
       
        ## Measured data       
        full_path = os.path.join(dossier, fichier)
        data = pd.read_csv(full_path, sep=',')
        f = data['Frequency']
        P = data['Power_Loss']

        ##Simulated data
        mfem_file = "Princeton_mfem"+mat+'_'+phi+".csv"
        path_sim = os.path.join(path, "data", "MFEM", mfem_file)
        data_sim = pd.read_csv(path_sim, sep=";")
        f_sim=data_sim['fc']
        P_sim = data_sim['P_tot']


        plt.figure()
        plt.plot(f,P, label="Measurements")
        plt.plot(f_sim, P_sim, label="MFEM")
        plt.grid()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power Loss (W/m3)')
        plt.title(f"{mat}, b_peak={phi} mT")
        plt.legend()

plt.show()

