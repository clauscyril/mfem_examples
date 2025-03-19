import matplotlib.pyplot as plt
import numpy as np


# Get mesh vertices coordinates
Z = []
liste = []
path = 'build/refined.mesh'
with open(path, 'r', encoding='utf-8') as f:
    liste = f.readlines()
    for i, e in enumerate(liste):
        if "vertices" in e :
            Z = liste[i+3:]
            break

Z = [float(e.replace("\n", "")) for e in Z]
# print(Z)


# Get solution
paths = ["build/sol_i.gf", "build/sol_r.gf"]
values = [[], []]

for k,path in enumerate(paths):
    with open(path, 'r', encoding='utf-8') as f:
        data = f.readlines()
        data = data[5:]
        for i in range(len(data)):
            values[k].append(float(data[i].replace("\n", '')))
                  
u_abs = []

for i in range(len(values[0])):
    # u_abs.append((values[0][i] ** 2 + values[1][i] ** 2)**0.5)
    u_abs.append(values[1][i])

n = len(u_abs)
e = 0.001

Z = np.array(Z)
u_abs = np.array(u_abs)


sorted_indices = np.argsort(Z)
Z = Z[sorted_indices] - e/2
u_abs = u_abs[sorted_indices]

f = 400
omega = 2*np.pi*f
mu = 4e-7*np.pi * 100
sigma = 10**7 

delta = np.sqrt(2/omega*mu*sigma)
A = np.sqrt(omega*sigma*mu*1j)




# u_ex = lambda z : 1/np.sinh(A*e) * (np.sinh(A*e/2 + A*z) + np.sinh(A*e/2 - A*z))
u_ex = lambda z : np.sinh(A*e/2)/np.sinh(A*e) * (np.exp(A*z) + np.exp(-A*z))
Z2 = np.linspace(-e/2, e/2, 100)
U_ex = np.array([np.real(u_ex(z)) for z in Z2])

# erreur = np.abs( (u_abs - np.real(u_ex(Z)))/np.real(u_ex(Z))) * 100
erreur = np.abs((u_abs - np.real(u_ex(Z))))

# print(Z)
# print(U_ex)

# print(delta*1000)


plt.figure()

plt.subplot(1,2,1)
plt.plot(1000*Z, u_abs, '--x', label='MFEM')
plt.plot(1000*Z2, U_ex, label="Solution exacte")
plt.xlabel("z (mm)")
plt.ylabel('H (A/m)')
plt.title("Résultats simulation MFEM")
plt.grid()
plt.legend()

plt.subplot(1,2,2)
plt.plot(1000*Z,erreur, '--x')
plt.title("Erreurs")
plt.xlabel("z(mm)")
plt.ylabel("H (A/m)")
plt.grid()

plt.show()


