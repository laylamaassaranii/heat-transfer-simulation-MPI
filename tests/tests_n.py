import math
import os

# =========================
# Common physical parameters
# =========================
Lx = 1.0
Ly = 1.0
alpha = 0.1
T_final = 0.5
dt_in = 0.0
bc_str = "NEUMANN"

# =========================
# Output folder
# =========================
data_dir = "data_tests_n"
os.makedirs(data_dir, exist_ok=True)

# =========================
# Helper: write one file
# =========================
def write_initial_file(Nx, Ny, filename):
    path = os.path.join(data_dir, filename)
    with open(path, "w") as f:
        # Header (must match C++ exactly)
        f.write(f"{Nx} {Ny}\n")
        f.write(f"{Lx} {Ly}\n")
        f.write(f"{alpha}\n")
        f.write(f"{T_final}\n")
        f.write(f"{dt_in}\n")
        f.write(f"{bc_str}\n")

        # Temperature field: Ny lines, each with Nx values (Gaussian hot spot)
        cx = (Nx - 1) / 2.0
        cy = (Ny - 1) / 2.0
        sigma = Nx / 10.0  # same “width” scaling as your example

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

    print(f"Wrote {path}  (Nx*Ny = {Nx*Ny}, values = {count})")


# =========================
# 1) Strong scaling tests
#    (same problem, more processes)
# =========================
strong_sizes = [
    (512, 512),
    (1024, 1024),
]

for Nx, Ny in strong_sizes:
    fname = f"initial_strong_{Nx}x{Ny}.txt"
    write_initial_file(Nx, Ny, fname)


# =========================
# 2) Weak scaling tests
#    (keep work per process ~ constant)
#    P = 1, 4, 9, 16
# =========================
weak_configs = [
    (1, 256, 256),    # P=1  -> ~256x256
    (4, 512, 512),    # P=4  -> ~512x512
    (9, 768, 768),    # P=9  -> ~768x768
    (16, 1024, 1024), # P=16 -> ~1024x1024
]

for P, Nx, Ny in weak_configs:
    fname = f"initial_weak_P{P}_{Nx}x{Ny}.txt"
    write_initial_file(Nx, Ny, fname)


# =========================
# 3) Optional: your original 200x200 test
# =========================
write_initial_file(200, 200, "initial_200x200_reference.txt")



"""
# strong scaling
time ./seq data_tests/initial_strong_512x512.txt seq_512
time mpiexec -n 4 ./par data_tests/initial_strong_512x512.txt par_512_np4

# weak scaling
time mpiexec -n 1  ./par data_tests/initial_weak_P1_256x256.txt  out_256_np1
time mpiexec -n 4  ./par data_tests/initial_weak_P4_512x512.txt  out_512_np4
...
"""