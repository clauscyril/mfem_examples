import matplotlib.pyplot as plt
import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))
path_file0 = os.path.join(path,"build", "data0.csv")
path_file1 = os.path.join(path,"build", "data1.csv")
path_file2 = os.path.join(path,"build", "data2.csv")
path_file3 = os.path.join(path,"build", "data3.csv")
path_file4 = os.path.join(path,"build", "data4.csv")


data0 = pd.read_csv(path_file0, sep=";")
data1 = pd.read_csv(path_file1, sep=";")
data2 = pd.read_csv(path_file2, sep=";")
data3 = pd.read_csv(path_file3, sep=";")
data4 = pd.read_csv(path_file4, sep=";")


Fc = data0['fc']/1000
Ploss = data0['Ploss']/1000

Fc1 = data1['fc']/1000
Ploss1 = data1['Ploss']/1000

Fc2 = data2['fc']/1000
Ploss2 = data2['Ploss']/1000

Fc3 = data3['fc']/1000
Ploss3 = data3['Ploss']/1000

Fc4 = data4['fc']/1000
Ploss4 = data4['Ploss']/1000

print(Ploss2)
# Création du graphique
plt.figure()
plt.plot(Fc, Ploss, '-*', label="Tau 0")
plt.plot(Fc1, Ploss1, label="Tau 1")
plt.plot(Fc2, Ploss2, '-*', label="Tau 2")
plt.plot(Fc3, Ploss3, label="Tau 3")
plt.plot(Fc4, Ploss4, label="Tau 4")

# Définition de l'échelle logarithmique pour l'axe y
plt.xscale('log')

# Ajout de labels et titre
plt.xlabel('Frequency (kHz)')
plt.ylabel('Losses (kW/m^3)')
plt.title('Ferrite N87')
plt.grid()
plt.legend()
# Affichage du graphique
plt.show()



# for i in range(4):
#     name = "data"
#     name += str(i) + "csv"
#     print(name)