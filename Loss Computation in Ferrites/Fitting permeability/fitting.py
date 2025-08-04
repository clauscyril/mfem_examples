import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Données d'exemple
w_data = np.linspace(-10, 10, 100)
true_a = 5
true_b = 2
y_data = true_a / np.sqrt(1 + (true_b**2) * w_data**2) + np.random.normal(0, 0.1, size=w_data.shape)

# Modèle
def f(w, a, b):
    return a / np.sqrt(1 + (b**2) * w**2)

# Ajustement
popt, pcov = curve_fit(f, w_data, y_data, p0=[1, 1])
fitted_a, fitted_b = popt

# Result
print(f"Ajustement: a = {fitted_a:.3f}, b = {fitted_b:.3f}")

# Visualisation
import matplotlib.pyplot as plt
plt.scatter(w_data, y_data, label='Données')
plt.plot(w_data, f(w_data, *popt), color='red', label='Ajustement')
plt.legend()
plt.show()
