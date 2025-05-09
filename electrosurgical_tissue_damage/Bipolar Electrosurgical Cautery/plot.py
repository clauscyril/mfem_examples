import matplotlib.pyplot as plt
import pandas as pd
import os


path_file = os.path.dirname(os.path.abspath(__file__))
path_file1 = os.path.join(path_file, "build", "data.csv")
path_file2 = os.path.join(path_file, "build", "data2.csv")

data = pd.read_csv(path_file1, sep=";")
data2 = pd.read_csv(path_file2, sep=";")


time = data["Time"]
temperature = data["Temperature"]
t_1mm = data["T_1mm"]

temperatur2 = data2['Temperature']

# plt.plot(time, temperature, label="Avec convexion")
# plt.plot(time, temperatur2, label="Sans convexion")

plt.plot(time, temperature, label="Au milieu")
plt.plot(time, t_1mm, label="À 1mm")

plt.grid()
plt.legend()
plt.show()