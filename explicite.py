import numpy as np
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
U = np.ones(shape = (Nt, Nz)) * 15

# Z axis (depth)
Z = np.linspace(0, L, Nz)

# T axis (time)
T = np.linspace(0, S, Nt)

# Initiale condition (temperature surface approximation)
def Uini(t):
    res = 15 - 10 * np.sin((2 * np.pi * t) / 12)
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
while err > 10**-4:

    maxiter -=1
    # Sauvegarde de la dernière temperature calculée
    Uold = U.copy()

    # Calcul de U au temps suivant pour chaque profondeur
    for t in range (1, Nt):
        U[t, 1: -1] = U[t-1, 1: -1] + (K * dt / (dz**2)) * (U[t-1, 0: -2] - 2 * U[t-1, 1: -1] + U[t-1, 2:])
    
    # Condition au bord (fond)
    U[:, -1] = U[:, -2]

    if maxiter == 0:
        print("Instable solution")
        raise SystemExit

    # Recherche de l'État Stationnaire Périodique
    err = np.linalg.norm(U - Uold)

    # Continuité entre le fin d'année précédente et le debut de la nouvelle année.
    U[0, :] = U[-1, :]


# End timer
end_time = time.time()

total_time = start_time - end_time

print("Total time execution: ", total_time)
