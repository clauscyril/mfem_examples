import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


path = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(path, "data", "15mT", "test2.csv")

data = pd.read_csv(path, sep=",")





B = data["B [mT]"] / 1000
H = data["H [A/m]"]
t = data["t"]


f = 433e3
T = 1/f
Ts = T/len(B)


height = 7.59e-3
w = 5.3e-3
b_peak = 10e-3
phi_peak = b_peak*w*height

# f = 100e3
# I_rms = 0.082/np.sqrt(2)
# nb_period = 1
# nb_iter = 1000 * nb_period
# Ts = nb_period/f/nb_iter

rho = 5.98e-2
sigma = 4.44e-1
eps = 2.48e-6
mu = 4300 * 4*np.pi*1e-7

tau = 1/(2*np.pi*1.8e6)


A1 = (rho * (sigma * Ts + 2*eps) + Ts)/(2*eps + Ts*sigma)
A2 = (rho * (sigma * Ts - 2*eps) + Ts)/(2*eps + Ts*sigma)
A3 = -(Ts*sigma -2*eps) / (2*eps + Ts*sigma)

B1 = 2*mu / (Ts + 2*tau)
B2 = -(2*mu) / (Ts + 2*tau)
B3 = -(Ts - 2*tau) / (Ts + 2*tau)

C1 = Ts*mu/(Ts+2*tau)
C2 = Ts*mu/(Ts+2*tau)
C3 = -(Ts-2*tau)/(Ts+2*tau)

D1 = Ts*mu/(Ts + tau)
D3 = tau/(Ts+tau)

new_B = []
Hnm1 = 0
Bnm1 = 0


for i in range(len(H)) :
    Hn = H[i]
    new_B.append(C1*Hn + C2* Hnm1 + C3*Bnm1)
    Hnm1 = Hn
    Bnm1 = new_B[i]

phi_nm1 = 0
new_H = []
for i in range(len(H)) :
    phi_n = B[i]
    new_H.append(1/D1 * (phi_n - D3 * phi_nm1))
    phiH_nm1 = new_H[i]
    phi_nm1 = phi_n; 

print(Ts)


plt.plot(t, B, label="B measured")
plt.plot(t, new_B, label="B computed from H")
plt.grid()
plt.legend()

plt.figure()
plt.plot(t, H, label="H measured")
plt.plot(t, new_H, label="H from B")
plt.legend()
plt.grid()

plt.show()