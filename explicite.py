import numpy as np
import matplotlib.pyplot as plt
import time

##  Initialisation of phyiscals parameters

# Total depth (meters)
L = 100

# Total time of one year (in seconds)
S = 365 * 24 * 60 * 60

# Dimension of the meshing
Nz = 400
Nt = 5000

# Initialisation of meshing (15 degres celsiuce everywhere)
U = np.ones(shape = (Nt + 1, Nz+ 1)) * 15

# Z axis (depth)
Z = np.linspace(0, L, Nz + 1)

# T axis (time)
T = np.linspace(0, S, Nt + 1)

# Time in months in order to create cycle 
Tmonths = np.array([ i+1 for i in range(0,12)])

# Initiale condition (temperature surface approximation)
def Uini(t):
    res = 15 - 10 * np.sin((2 * np.pi * t) / S)
    return res

# Finale condition 
def Ufin(t,z):
    res = 0
    return res

# Pas d'espace 
dz = L / Nz

# Pas de temps 
dt = S / Nt

# Coeff of difusivity
K = 10**-6

## Calcul par differences finies

# Timer
start_time = time.time()

# Calcul de U à la surface pour tout les t
U[:, 0] = Uini(T)

# Condition d'arret si non convergence
maxiter = 500

# On met a jour chaque année la température initial avec la temperature final de l'an précédent
# Obtenir un 
err = 1
y = 0
R = []
Y = []

# Doit etre inférieur à 1/2 sinon instabilités
r = (K * dt / (dz**2))
r = 0.499
while err > 10**-4:

    maxiter -=1
    # Sauvegarde de la dernière temperature calculée
    Uold = U.copy()

    # Calcul de U au temps suivant pour chaque profondeur
    for t in range (1, Nt):
        U[t, 1: -1] = U[t-1, 1: -1] + r * (U[t-1, 0: -2] - 2 * U[t-1, 1: -1] + U[t-1, 2:])
    
    # Condition au bord (fond)
    U[:, -1] = U[:, -2]

    if maxiter == 0:
        print("Instable solution")
        raise SystemExit

    # Continuité entre le fin d'année précédente et le debut de la nouvelle année.
    U[0, :] = U[-1, :]

    # Recherche de l'État Stationnaire Périodique
    err = np.linalg.norm(U - Uold)
    R.append(err)
    y += 1
    Y.append(y)
   
# End timer
end_time = time.time()
total_time = start_time - end_time

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

x = 2 
y2 = 3
ax1.plot(Y, R, '-', linewidth=2)
# S/12 = nbr seconde in a months = 2628000xlim=(0, 5), ylim=(-5, 35 )
for d in range(0,100,15):
    ax2.plot(T /  (S / 12), U[:, d],'-', label = f"{d} meters", linewidth=2)
ax3.plot(x, y, '-', linewidth=2.0)
ax4.plot(x, y2 - 2.5, 'o-', linewidth=2)

ax1.set_title("error evolution", color='red', fontsize=10)
ax1.set(xlim=(0, 2), xlabel = "time in year", ylabel = "error")
ax2.set_title("1D temperature variation", color='red', fontsize=10)
ax2.set(xlabel = "time (months)", ylabel = "temperature C°")
ax2.legend()
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()

img = ax3.imshow(U.transpose(), 
                 aspect='auto', 
                 extent=[0, 12, 100, 0], 
                 cmap='jet',        # 'jet' est la couleur classique bleu-rouge (comme Matlab)
                 interpolation='nearest')


ax3.set_title("Variation de Temperature (imagesc - imshow)")
ax3.set_xlabel("Temps (Mois)")
ax3.set_ylabel("Profondeur (m)")
plt.colorbar(img, ax=ax3, label="Température (°C)")
X_grid, Y_grid = np.meshgrid(T / (S / 12), Z)

# On dessine les contours remplis
# 20 = nombre de niveaux de couleurs (plus il est grand, plus c'est lisse)
cnt = ax4.contourf(X_grid, Y_grid, U.T, 20, cmap='jet')

# Important : Inverser l'axe Y pour avoir la surface (0m) en HAUT
ax4.invert_yaxis()

ax4.set_title("Variation de Temperature (contourf)")
ax4.set_xlabel("Temps (Mois)")
ax4.set_ylabel("Profondeur (m)")
plt.colorbar(cnt, ax=ax4, label="Température (°C)")

plt.tight_layout()

plt.show()

A = 1