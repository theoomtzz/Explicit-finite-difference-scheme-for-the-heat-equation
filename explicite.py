import numpy as np

##  Initialisation of phyiscals parameters

# Total depth (meters)
L = 100

# Total time of one year (in seconds)
S = 365 * 24 * 60 * 60

# Dimension of the meshing
Nz = 400
Nt = 5000

# Initialisation of meshing
data = np.zeros(shape = (Nt, Nz))

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




