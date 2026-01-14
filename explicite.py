import numpy as np
import matplotlib.pyplot as plt
import time

# ==============================================================================
# PARAMETERS & CONFIGURATION
# ==============================================================================

# Physical constants
# ------------------
# L : Total depth of the spatial domain in meters
L = 100

# S : One year expressed in seconds (Temporal period)
S = 365 * 24 * 60 * 60

# K : Thermal Diffusivity Constant (m^2/s)
K = 10**-6

# Discretization
# --------------
# Nz : Number of spatial steps (Spatial Mesh)
# Nt : Number of time steps (Temporal Mesh)
Nz = 400
Nt = 5000

# Grid Initialization
# -------------------
# Z : Spatial Grid Vector (meters)
Z = np.linspace(0, L, Nz + 1)

# T : Temporal Grid Vector (seconds)
T = np.linspace(0, S, Nt + 1)

# U : State Matrix (Temperature field)
# Initialized at 15°C (Flat start)
U = np.ones(shape=(Nt + 1, Nz + 1)) * 15

# ==============================================================================
# BOUNDARY CONDITIONS
# ==============================================================================

def Uini(t):
    """
    Computes the Dirichlet boundary condition at the surface (z=0).
    Simulates seasonal temperature variation.

    Parameters:
    -----------
    t : Current simulation time in seconds (float or numpy array)

    Returns:
    --------
    temp : Surface temperature in Celsius (float or numpy array)
    """
    # Sinusoidal fluctuation around 15°C with amplitude of 10°C
    res = 15 - 10 * np.sin((2 * np.pi * t) / S)
    return res

def Ufin(t, z):
    """
    Placeholder for Bottom Boundary Condition.
    (Not strictly used in the explicit loop below as we force Neumann BC manually).
    
    Parameters:
    -----------
    t : Time (float)
    z : Depth (float)

    Returns:
    --------
    val : Boundary value (float)
    """
    return 0

# ==============================================================================
# NUMERICAL SOLVER (EXPLICIT SCHEME)
# ==============================================================================

# Pre-compute mesh steps
# ----------------------
dz = L / Nz  # Spatial step size
dt = S / Nt  # Time step size

# Stability Factor (CFL Condition)
# ------------------------------
# r must be < 0.5 for the explicit scheme to be stable.
# We hardcode it to 0.499 to maximize step size while maintaining stability.
r = (K * dt / (dz**2))
r = 0.499 

# Convergence Loop Initialization
# -------------------------------
start_time = time.time()

# Apply Surface Boundary Condition for the entire time vector
U[:, 0] = Uini(T)

maxiter = 500  # Max years to simulate
err = 1        # Error initialization
y = 0          # Year counter
R = []         # Error history
Y = []         # Iteration history

while err > 10**-4:
    maxiter -= 1
    
    # Store previous state to compute convergence metric later
    Uold = U.copy()

    # Time-Stepping Loop
    # ------------------
    # Vectorized Forward Euler / Central Difference in Space
    # u_i^{n+1} = u_i^n + r * (u_{i-1}^n - 2u_i^n + u_{i+1}^n)
    for t in range(1, Nt):
        U[t, 1:-1] = U[t-1, 1:-1] + r * (U[t-1, 0:-2] - 2 * U[t-1, 1:-1] + U[t-1, 2:])
    
    # Boundary Condition at Bottom (z=L)
    # ----------------------------------
    # Neumann Condition: Adiabatic / No-flux (dT/dz = 0) implies U[-1] = U[-2]
    U[:, -1] = U[:, -2]

    # Stability Check
    if maxiter == 0:
        print("Error: Max iterations reached. Solution might be unstable.")
        raise SystemExit

    # Continuity Constraint
    # ---------------------
    # The end state of the current year becomes the start state of the next year
    U[0, :] = U[-1, :]

    # Convergence Check
    # -----------------
    # Compute L2 Norm of the difference between two consecutive years
    err = np.linalg.norm(U - Uold)
    R.append(err)
    y += 1
    Y.append(y)
   
end_time = time.time()

# ==============================================================================
# VISUALIZATION
# ==============================================================================

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# 1. Convergence Plot
# -------------------
ax1.plot(Y, R, '-', linewidth=2)
ax1.set_title("Error Evolution (Convergence)", color='red', fontsize=10)
ax1.set(xlim=(0, len(Y)), xlabel="Year Iteration", ylabel="Error Norm")
ax1.grid(True, linestyle='--', alpha=0.6)

# 2. Temperature Profiles (1D Slices)
# -----------------------------------
# Normalize time to months for readability
T_months = T / (S / 12)

for d in range(0, 100, 15):
    # Plotting temperature vs time for specific depths
    ax2.plot(T_months, U[:, d], '-', label=f"{d} meters", linewidth=2)

ax2.set_title("1D Temperature Variation", color='red', fontsize=10)
ax2.set(xlabel="Time (months)", ylabel="Temperature (°C)")
ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
ax2.grid(True, linestyle='--', alpha=0.6)

# 3. Heatmap (Raw Data)
# ---------------------
# Transpose U to align axes: Y=Depth, X=Time
img = ax3.imshow(U.transpose(), 
                 aspect='auto', 
                 extent=[0, 12, 100, 0], 
                 cmap='jet',        
                 interpolation='nearest')

ax3.set_title("Temperature Heatmap (imshow)")
ax3.set_xlabel("Time (Months)")
ax3.set_ylabel("Depth (m)")
plt.colorbar(img, ax=ax3, label="Temperature (°C)")

# 4. Contour Plot (Interpolated)
# ------------------------------
X_grid, Y_grid = np.meshgrid(T_months, Z)
cnt = ax4.contourf(X_grid, Y_grid, U.T, 20, cmap='jet')
ax4.invert_yaxis() # Ensure depth 0 is at the top

ax4.set_title("Temperature Isocontours (contourf)")
ax4.set_xlabel("Time (Months)")
ax4.set_ylabel("Depth (m)")
plt.colorbar(cnt, ax=ax4, label="Temperature (°C)")

plt.tight_layout()
plt.show()