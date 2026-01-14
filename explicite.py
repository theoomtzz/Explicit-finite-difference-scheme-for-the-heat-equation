import numpy as np
import matplotlib.pyplot as plt
import time

# ==============================================================================
# PHYSICAL PARAMETERS & DOMAIN DISCRETIZATION
# ==============================================================================

# Total depth of the spatial domain (L) in meters
L = 100

# Temporal Period: One year expressed in seconds (T_cycle)
S = 365 * 24 * 60 * 60

# Grid Resolution:
# Nz: Number of spatial steps (Spatial Mesh)
# Nt: Number of time steps (Temporal Mesh)
Nz = 400
Nt = 5000

# State Matrix Initialization:
# Field U(t,z) initialized at 15°C (Flat start)
U = np.ones(shape = (Nt + 1, Nz+ 1)) * 15

# Spatial Grid Vector (Z axis)
Z = np.linspace(0, L, Nz + 1)

# Temporal Grid Vector (T axis)
T = np.linspace(0, S, Nt + 1)

# Time mapping for plotting (Months) - Auxiliary array for cycling
Tmonths = np.array([ i+1 for i in range(0,12)])

# ------------------------------------------------------------------------------
# BOUNDARY CONDITIONS FUNCTIONS
# ------------------------------------------------------------------------------

# Dirichlet Boundary Condition at Surface (z=0)
# Simulates seasonal temperature variation
def Uini(t):
    # Initial condition (temperature surface approximation)
    res = 15 - 10 * np.sin((2 * np.pi * t) / S)
    return res

# Placeholder for Bottom Boundary Condition (Not used in explicit loop logic below)
def Ufin(t,z):
    res = 0
    return res

# ------------------------------------------------------------------------------
# NUMERICAL SCHEME PARAMETERS
# ------------------------------------------------------------------------------

# Spatial Step size (Delta z) - Pas d'espace 
dz = L / Nz

# Time Step size (Delta t) - Pas de temps 
dt = S / Nt

# Thermal Diffusivity Constant (K) - Coeff of difusivity
K = 10**-6

# ==============================================================================
# FINITE DIFFERENCE SOLVER (EXPLICIT SCHEME)
# ==============================================================================

# Performance Timer Start
start_time = time.time()

# Apply Surface Boundary Condition for all time steps (Vectorized)
# Calcul de U à la surface pour tout les t
U[:, 0] = Uini(T)

# Convergence criteria (Max Iterations) - Condition d'arret
maxiter = 500

# Convergence Loop Setup
# We iterate year-by-year until a Periodic Steady State is reached.
# On met a jour chaque année la température initial avec la temperature final de l'an précédent
err = 1
y = 0
R = [] # Error history container
Y = [] # Iteration history container

# STABILITY CONDITION (CFL - Courant-Friedrichs-Lewy)
# Must be < 0.5 for explicit scheme stability.
# Hardcoded to 0.499 to maximize step size while remaining stable.
# Doit etre inférieur à 1/2 sinon instabilités
r = (K * dt / (dz**2))
r = 0.499 

# Main Solver Loop
while err > 10**-4:

    maxiter -=1
    
    # Snapshot of the previous state for convergence check
    # Sauvegarde de la dernière temperature calculée
    Uold = U.copy()

    # Time-Stepping Loop (Forward Euler / Central Difference in Space)
    # Calcul de U au temps suivant pour chaque profondeur
    for t in range (1, Nt):
        # Explicit update rule: u_i^{n+1} = u_i^n + r * (u_{i-1}^n - 2u_i^n + u_{i+1}^n)
        U[t, 1: -1] = U[t-1, 1: -1] + r * (U[t-1, 0: -2] - 2 * U[t-1, 1: -1] + U[t-1, 2:])
    
    # Neumann Boundary Condition at Bottom (z=L)
    # Adiabatic / No-flux condition: dT/dz = 0 => U[-1] = U[-2]
    # Condition au bord (fond)
    U[:, -1] = U[:, -2]

    # Safety break for non-convergence
    if maxiter == 0:
        print("Instable solution / Max iterations reached")
        raise SystemExit

    # Continuity Constraint:
    # The end state of the current year becomes the start state of the next iteration
    # Continuité entre le fin d'année précédente et le debut de la nouvelle année.
    U[0, :] = U[-1, :]

    # Check for Periodic Steady State (L2 Norm of difference)
    # Recherche de l'État Stationnaire Périodique
    err = np.linalg.norm(U - Uold)
    R.append(err)
    y += 1
    Y.append(y)
   
# Performance Timer End
end_time = time.time()
total_time = start_time - end_time

# ==============================================================================
# DATA VISUALIZATION / REPORTING
# ==============================================================================

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

# Placeholder data for unused subplots (debug traces)
x = 2 
y2 = 3

# Plot 1: Convergence History
ax1.plot(Y, R, '-', linewidth=2)
ax1.set_title("Error Evolution (Convergence)", color='red', fontsize=10)
ax1.set(xlim=(0, 2), xlabel = "Year Iteration", ylabel = "Error Norm")

# Plot 2: Temperature Profiles at specific depths
# Normalizing time to Months for X-axis (S/12 = seconds in a month)
for d in range(0,100,15):
    ax2.plot(T / (S / 12), U[:, d],'-', label = f"{d} meters", linewidth=2)

ax2.set_title("1D Temperature Variation", color='red', fontsize=10)
ax2.set(xlabel = "Time (months)", ylabel = "Temperature (°C)")
ax2.legend()
# Moves legend outside the plot area for clarity
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Plot 3 & 4 Placeholders (Leftover from debugging)
ax3.plot(x, y, '-', linewidth=2.0)
ax4.plot(x, y2 - 2.5, 'o-', linewidth=2)

plt.tight_layout()

# --- Heatmap Visualization (2D Field) ---

# Imshow: Raw matrix visualization
# Transpose U to map Depth to Y-axis and Time to X-axis
img = ax3.imshow(U.transpose(), 
                 aspect='auto', 
                 extent=[0, 12, 100, 0], 
                 cmap='jet',        
                 interpolation='nearest')

ax3.set_title("Variation de Temperature (imagesc - imshow)")
ax3.set_xlabel("Temps (Mois)")
ax3.set_ylabel("Profondeur (m)")
plt.colorbar(img, ax=ax3, label="Température (°C)")

# --- Contour Visualization (Interpolated Field) ---

# Meshgrid generation for contour plotting
X_grid, Y_grid = np.meshgrid(T / (S / 12), Z)

# Contourf: Filled contour plot for smoother gradient visualization
cnt = ax4.contourf(X_grid, Y_grid, U.T, 20, cmap='jet')

# Invert Y-axis to represent depth correctly (0m at surface/top)
ax4.invert_yaxis()

ax4.set_title("Variation de Temperature (contourf)")
ax4.set_xlabel("Temps (Mois)")
ax4.set_ylabel("Profondeur (m)")
plt.colorbar(cnt, ax=ax4, label="Température (°C)")

plt.tight_layout()

plt.show()

A = 1