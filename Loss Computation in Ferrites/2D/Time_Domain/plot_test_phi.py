import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

f = 1000e3
mu = 4300.0 * 4e-7 * np.pi
omega = 2*np.pi*f;           
tau = 1./(2.*np.pi*1.8e6);   
mu_eq = mu / (1. + tau * 1j * omega)
phase = np.abs(mu_eq)

X = np.linspace(0,5/f, 1000)
phi_phased = np.sin(2*np.pi*f*X - phase)

path = r"C:\Users\cyril\Projets\mfem_examples\Loss Computation in Ferrites\2D\Time_Domain\build\test.csv"
data = pd.read_csv(path, sep=";")
t = data['t']
phi = data['phi']
phi_H = data['phiH']

plt.plot(t,phi_H)
plt.grid()
plt.show()