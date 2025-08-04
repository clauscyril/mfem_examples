import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "data", "15mT")

list_f = [100, 150, 200, 250, 300, 351]
abs_mu = []


for f in list_f:
    data_path = os.path.join(path, str(f) + "kHz.csv")
    data = pd.read_csv(data_path, sep=",")

    B = data["B [mT]"]
    H = data["H [A/m]"]

    ratio = max(B)/max(H) * 1e-3/(4e-7*np.pi)
    print(ratio)
    abs_mu.append(ratio)

plt.plot(list_f, abs_mu)
plt.show()