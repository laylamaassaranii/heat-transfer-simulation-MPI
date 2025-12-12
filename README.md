# Heat Transfer Simulation - MPI vs Sequential

A 2D heat transfer simulation comparing parallel (MPI) and sequential execution methods. This project implements the 2D heat equation using finite difference methods to model temperature distribution over time.

## Contributors

- Layla Maassarani (201865)
- Ali Joumaa (211058)
- Neamat Khatib (220966)
- Farouk Tannir (201153)

## Project Overview

This project simulates heat transfer in a 2D grid using the heat equation:

$$\frac{\partial T}{\partial t} = \alpha \nabla^2 T$$

where:

- $T$ is the temperature field
- $\alpha$ is the thermal diffusivity
- $\nabla^2$ is the Laplacian operator

The simulation is implemented in two ways:

1. **Sequential** - Single-threaded CPU execution
2. **Parallel** - Distributed computing using MPI (Message Passing Interface)

## Project Structure

```
heat-transfer-simulation-MPI/
├── sequential/
│   ├── sequential.cpp           # Sequential implementation
│   ├── initial_conditions.txt   # Input for sequential version
│   └── heat_output.csv          # Output from sequential run
├── parallele/
│   ├── parallel.cpp             # MPI parallel implementation
│   ├── initial_conditions_large.txt  # Input for parallel version
│   └── heat_mpi_large_np4_2_final.csv # Output from parallel run
├── plots/                       # Output directory for visualizations
├── generate_ic_large.py         # Script to generate large initial conditions
├── visualization.py             # Script to visualize simulation results
└── README.md                    # This file
```

## Implementation Details

### Sequential Version ([sequential/sequential.cpp](sequential/sequential.cpp))

- **Method**: Standard finite difference method with explicit time-stepping
- **Boundary Conditions**: Dirichlet boundary conditions (fixed temperatures at edges)
- **Grid Update**: Updates entire grid at each time step sequentially
- **Best for**: Small to medium grid sizes (e.g., 50×50 to 100×100)

### Parallel Version ([parallele/parallel.cpp](parallele/parallel.cpp))

- **Method**: Domain decomposition using MPI
- **Parallelization Strategy**:
  - 2D Cartesian topology for process distribution
  - Each MPI process handles a subdomain of the grid
  - Ghost cell exchange between neighboring processes
- **Boundary Conditions**: Supports both Dirichlet and Neumann boundary conditions
- **Communication**: Non-blocking sends/receives for ghost cell updates
- **Best for**: Large grid sizes (e.g., 200×200 and above)

**Key Features:**

- Dynamic domain decomposition using `MPI_Dims_create()`
- Halo/ghost cells for boundary communication
- Stable time-step calculation based on CFL condition
- Efficient gather operation for final output

## Input File Format

Both implementations read initial conditions from text files with the following format:

```
Nx Ny                    # Grid dimensions
Lx Ly                    # Physical domain size
alpha                    # Thermal diffusivity
T_final                  # Final simulation time
dt                       # Time step (0.0 for auto-calculation)
BOUNDARY_TYPE            # DIRICHLET or NEUMANN (parallel only)
[Temperature values]     # Ny lines × Nx values per line
```

**Example:**

```
50 50
1.0 1.0
0.1
0.5
0.0
DIRICHLET
20.0 20.0 20.0 ... (temperature field)
```

## Usage

### 1. Generate Initial Conditions

For large-scale simulations, generate initial conditions using:

```bash
python generate_ic_large.py
```

This creates a `initial_conditions_large.txt` file with a Gaussian temperature distribution:

- Default: 200×200 grid
- Hot spot at center with Gaussian decay
- Parameters: Lx=Ly=1.0, α=0.1, T_final=0.5

### 2. Compile the Programs

**Sequential:**

```bash
g++ sequential/sequential.cpp -o sequential.exe -O3
```

**Parallel (requires MPI):**

```bash
mpic++ parallele/parallel.cpp -o parallel.exe -O3
```

### 3. Run Simulations

**Sequential:**

```bash
./sequential.exe sequential/initial_conditions.txt sequential/heat_output
```

**Parallel (with 4 processes):**

```bash
mpirun -np 4 ./parallel.exe parallele/initial_conditions_large.txt parallele/heat_mpi_large_np4_2
```

The output will be saved as `<output_prefix>_final.csv` for parallel or `<output_prefix>.csv` for sequential.

### 4. Visualize Results

Generate heatmap visualizations of the final temperature field:

```bash
python visualization.py sequential/heat_output.csv plots
python visualization.py parallele/heat_mpi_large_np4_2_final.csv plots
```

The visualization script:

- Reads CSV output files
- Creates a heatmap using matplotlib
- Saves figures with timestamps to the `plots/` directory
- Uses the "viridis" colormap for temperature visualization

## Algorithm Details

### Finite Difference Discretization

The 2D heat equation is discretized using central differences:

$$T_{i,j}^{n+1} = T_{i,j}^n + \alpha \Delta t \left[\frac{T_{i+1,j}^n - 2T_{i,j}^n + T_{i-1,j}^n}{\Delta x^2} + \frac{T_{i,j+1}^n - 2T_{i,j}^n + T_{i,j-1}^n}{\Delta y^2}\right]$$

### Stability Condition

The time step must satisfy the CFL (Courant-Friedrichs-Lewy) condition for stability:

$$\Delta t \leq \frac{1}{2\alpha\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right)}$$

The code automatically adjusts the time step to 90% of this maximum value if needed.

### Parallel Domain Decomposition

The parallel implementation:

1. Divides the grid into subdomains (one per MPI process)
2. Each process computes its local subdomain
3. Ghost cells are exchanged with neighbors at each time step
4. Process 0 gathers all results and writes the final output

## Performance Considerations

**Sequential Version:**

- Simple, easy to debug
- No communication overhead
- Memory efficient for small problems
- O(Nx × Ny × Nt) computational complexity

**Parallel Version:**

- Scales to large problems
- Communication overhead for ghost cell exchange
- Speedup depends on grid size and number of processes
- Optimal for problems where computation >> communication

## Dependencies

**For C++ code:**

- C++11 or later (tested with C++17)
- MPI library (OpenMPI, MPICH, or Intel MPI) for parallel version

**For Python scripts:**

- Python 3.6+
- NumPy
- Matplotlib

Install Python dependencies:

```bash
pip install numpy matplotlib
```

## Example Workflow

Complete workflow for running and visualizing a simulation:

```bash
# Generate large initial conditions
python generate_ic_large.py

# Compile both versions
g++ sequential/sequential.cpp -o sequential.exe -O3
mpic++ parallele/parallel.cpp -o parallel.exe -O3

# Run sequential version
./sequential.exe sequential/initial_conditions.txt sequential/heat_output

# Run parallel version with 4 processes
mpirun -np 4 ./parallel.exe parallele/initial_conditions_large.txt parallele/heat_mpi_large_np4_2

# Visualize both results
python visualization.py sequential/heat_output.csv plots
python visualization.py parallele/heat_mpi_large_np4_2_final.csv plots
```

## Output Files

- **CSV Files**: Comma-separated temperature values (Ny rows × Nx columns)
- **Plot Images**: PNG files with timestamp in `plots/` directory
- **Format**: Each row represents a horizontal line of the grid (j-index), each column is a vertical position (i-index)

## Notes

- The parallel version supports both Dirichlet (fixed) and Neumann (zero-gradient) boundary conditions
- The sequential version currently only supports Dirichlet boundary conditions
- Initial conditions can be customized by modifying `generate_ic_large.py`
- For best parallel performance, use grid sizes that are evenly divisible by the number of processes
- The visualization colormap can be changed in `visualization.py` (currently uses "viridis")
