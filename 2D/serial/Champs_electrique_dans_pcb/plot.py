import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
import os

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, "data", "sol.csv")

data = pd.read_csv(data_path)
# print(data)
x, y, potential = data["x"], data["y"], data["potentiel"]
triang = tri.Triangulation(x,y)

plt.figure(figsize=(8, 6))
contour = plt.tricontour(triang, potential, levels=30, cmap="coolwarm")
plt.colorbar(label="Potentiel")
plt.triplot(triang, 'k-', alpha=0.1)  # Tracer le maillage en noir
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Lignes d'équipotentielles")
plt.show()

