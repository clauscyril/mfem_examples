import matplotlib.pyplot as plt
import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))
path_file0_2D = os.path.join(path,"2D", "build", "data0.csv")
path_file1_2D = os.path.join(path,"2D","build", "data1.csv")
path_file2_2D = os.path.join(path,"2D","build", "data2.csv")
path_file3_2D = os.path.join(path,"2D","build", "data3.csv")
path_file4_2D = os.path.join(path,"2D","build", "data4.csv")


data0_2D = pd.read_csv(path_file0_2D, sep=";")
data1_2D = pd.read_csv(path_file1_2D, sep=";")
data2_2D= pd.read_csv(path_file2_2D, sep=";")
data3_2D = pd.read_csv(path_file3_2D, sep=";")
data4_2D = pd.read_csv(path_file4_2D, sep=";")


Fc_2D = data0_2D['fc']/1000
Ploss_2D = data0_2D['Ploss']/1000

Fc1_2D = data1_2D['fc']/1000
Ploss1_2D = data1_2D['Ploss']/1000

Fc2_2D = data2_2D['fc']/1000
Ploss2_2D = data2_2D['Ploss']/1000

Fc3_2D = data3_2D['fc']/1000
Ploss3_2D = data3_2D['Ploss']/1000

Fc4_2D = data4_2D['fc']/1000
Ploss4_2D = data4_2D['Ploss']/1000


path_file0_3D = os.path.join(path,"3D", "build", "data0.csv")
path_file1_3D = os.path.join(path,"3D","build", "data1.csv")
path_file2_3D = os.path.join(path,"3D","build", "data2.csv")
path_file3_3D = os.path.join(path,"3D","build", "data3.csv")
path_file4_3D = os.path.join(path,"3D","build", "data4.csv")


data0_3D = pd.read_csv(path_file0_3D, sep=";")
data1_3D = pd.read_csv(path_file1_3D, sep=";")
data2_3D= pd.read_csv(path_file2_3D, sep=";")
data3_3D = pd.read_csv(path_file3_3D, sep=";")
data4_3D = pd.read_csv(path_file4_3D, sep=";")


Fc_3D = data0_3D['fc']/1000
Ploss_3D = data0_3D['Ploss']/1000

Fc1_3D = data1_3D['fc']/1000
Ploss1_3D = data1_3D['Ploss']/1000

Fc2_3D = data2_3D['fc']/1000
Ploss2_3D = data2_3D['Ploss']/1000

Fc3_3D = data3_3D['fc']/1000
Ploss3_3D = data3_3D['Ploss']/1000

Fc4_3D = data4_3D['fc']/1000
Ploss4_3D = data4_3D['Ploss']/1000



# Création du graphique
plt.figure()
# plt.plot(Fc_2D, Ploss_2D, '--*', label="Tau 0")
# plt.plot(Fc1_2D, Ploss1_2D, '--*', label="Tau 1")
# plt.plot(Fc2_2D, Ploss2_2D, '--*', label="Tau 2")
# plt.plot(Fc3_2D, Ploss3_2D, '--*', label="Tau 3")
# plt.plot(Fc4_2D, Ploss4_2D, '--*', label="Tau 4")


plt.plot(Fc4_2D, Ploss4_2D, '--*', label="fc_mu = 10e6 2D")
plt.plot(Fc4_3D, Ploss4_3D, '--*', label="fc_mu = 10e6 3D")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')
# plt.yscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.title('Ferrite N87')
plt.grid()
plt.legend()
# Affichage du graphique
plt.show()
