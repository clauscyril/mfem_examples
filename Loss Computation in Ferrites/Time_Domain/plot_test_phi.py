import numpy as np
import matplotlib.pyplot as plt
import pandas as pd





height = 7.59e-3
w = 5.3e-3
b_peak = 10e-3
phi_peak = b_peak*w*height

f = 100e3
I_rms = 0.082/np.sqrt(2)
nb_period = 1
nb_iter = 1000 * nb_period
Ts = nb_period/f/nb_iter

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

def phi_sine_func(t) :
    return phi_peak *np.sin(2 * np.pi * f * t)


PhiH = [0]
Phi = [0]
Phi2 = [0]
PhiH2 = [0]
phi_n = 0
phiH_n = 0
phiH_nm1 = 0
phi_nm1 = 0

print("C1 = ", C1)
print("C2 = ", C2)
print("C3 = ", C3)


t = 0
for i in range(nb_iter):
    t+= Ts
    phi_n = phi_sine_func(t)
    phiH_n = 1/C1*(phi_n - C2 * phiH_nm1 - C3*phi_nm1)
    # print(f"phiH_n = 1//{C1} * ({phi_n} - {C2}*{phiH_nm1} - {C3}*{phi_nm1})")
    Phi2.append(C1*phiH_n + C2*phiH_nm1 + C3*Phi2[-1])
    PhiH2.append(1/D1 * (phi_n - D3 * phi_nm1))
    phiH_nm1 = phiH_n
    phi_nm1 = phi_n 
    PhiH.append(phiH_n)
    Phi.append(phi_n)

# t = 0
# for i in range(nb_iter):
#     t+= Ts
#     phiH_n = phi_sine_func(t)
#     phi_n = C1*phiH_n + C2*phiH_nm1 + C3*phi_nm1
#     phiH_nm1 = phiH_n
#     phi_nm1 = phi_n 
#     PhiH.append(phiH_n)
#     Phi.append(phi_n)



T = np.linspace(0,nb_iter*Ts,nb_iter + 1)
print(len(T))


list_f = np.linspace(10e3,2e6, 10000)
mu_complex = mu/(1+tau*2*np.pi*1j*list_f)


plt.figure()
plt.plot(T, Phi, label="Flux de B")
# plt.plot(T, Phi2, '--', label="Flux de B")
plt.legend()
plt.figure()
plt.plot(T, PhiH, label="Flux de H")
plt.plot(T, PhiH2, label="test ")
plt.legend()
plt.figure()
plt.plot(list_f, np.abs(mu_complex), label="Flux de H")
plt.loglog()
plt.grid()
plt.legend()
plt.show()
