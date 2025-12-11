import math
import os

# --- PARAMETERS ---
Nx = 200     # you can increase later to 800, 1000, etc.
Ny = 200

Lx = 1.0
Ly = 1.0
alpha = 0.1
T_final = 0.5
dt_in = 0.0        # let C++ choose stable dt
bc_str = "DIRICHLET"  # or "NEUMANN"

output_file = "initial_conditions_large.txt"

with open(output_file, "w") as f:
    # Header: must match your C++ exactly
    f.write(f"{Nx} {Ny}\n")
    f.write(f"{Lx} {Ly}\n")
    f.write(f"{alpha}\n")
    f.write(f"{T_final}\n")
    f.write(f"{dt_in}\n")
    f.write(f"{bc_str}\n")

    # Temperature field: Ny lines, each with Nx values
    cx = (Nx - 1) / 2.0
    cy = (Ny - 1) / 2.0
    sigma = Nx / 10.0

    count = 0
    for j in range(Ny):
        row_vals = []
        for i in range(Nx):
            dx = i - cx
            dy = j - cy
            r2 = dx * dx + dy * dy
            T = 20.0 + 80.0 * math.exp(-r2 / (2.0 * sigma * sigma))
            row_vals.append(f"{T:.6f}")
            count += 1
        f.write(" ".join(row_vals) + "\n")

print(f"Wrote {output_file} in {os.getcwd()}")
print(f"Nx * Ny = {Nx * Ny}, values written = {count}")
