import numpy as np
import matplotlib.pyplot as plt


f = np.linspace(100e3,2e6, 50)
rho = 5.98e-2
sigma = 4.44e-1
eps = 2.48e-6
mu = 4300.0 * 4e-7 * np.pi

rho_eq =  rho + 1. /(sigma + 1j * 2*np.pi*f * eps); 
mu_eq = mu / (1. + 1/(2*np.pi*1.8e6) * 1j * np.pi*f)

plt.figure()
plt.plot(f/1000, np.real(rho_eq))
plt.plot(f/1000, np.real(rho_eq) - np.imag(rho_eq))
plt.figure()
plt.plot(f/1000, np.imag(rho_eq))
plt.plot(f/1000, np.real(rho_eq) + np.imag(rho_eq))
plt.figure()
plt.plot(f/1000, np.real(mu_eq))
plt.plot(f/1000, np.real(mu_eq) - np.imag(mu_eq))
plt.figure()
plt.plot(f/1000, np.imag(mu_eq))
plt.plot(f/1000, np.real(mu_eq) + np.imag(mu_eq))




plt.show()
