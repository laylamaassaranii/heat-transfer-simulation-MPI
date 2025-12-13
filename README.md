# Heat Transfer Simulation - MPI vs Sequential

A 2D heat transfer simulation comparing parallel (MPI) and sequential execution methods. This project implements the 2D heat equation using finite difference methods to model temperature distribution over time.

## Contributors

- Layla Maassarani (201865)
- Ali Joumaa (211058)
- Neamat Khatib (220966)
- Farouk Tannir (201153)

## Documentation

[![PDF Report](https://img.shields.io/badge/Report-PDF-red?style=for-the-badge&logo=adobe)](./2d-heat-transfer-report.pdf)

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
│   ├── initial_conditions.txt   # Sample input file
│   └── heat_output.csv          # Sample output file
├── parallele/
│   ├── parallel.cpp             # MPI parallel implementation
│   ├── initial_conditions_large.txt  # Sample input file
│   └── heat_mpi_large_np4_2_final.csv # Sample output file
├── csv/
│   ├── seq/                     # Output directory for sequential runs
│   │   └── *.csv                # Sequential simulation results
│   └── par/                     # Output directory for parallel runs
│       └── *.csv                # Parallel simulation results
├── initial_data/
│   ├── data_tests_d/            # Dirichlet test initial conditions
│   └── data_tests_n/            # Neumann test initial conditions
├── tests/
│   ├── tests_d.py               # Generate Dirichlet test cases
│   └── tests_n.py               # Generate Neumann test cases
├── plots/                       # Output directory for visualizations
├── generate_ic_large.py         # Script to generate large initial conditions
├── visualization.py             # Script to visualize simulation results
└── README.md                    # This file
```

## Implementation Details

### Sequential Version ([sequential/sequential.cpp](sequential/sequential.cpp))

- **Method**: Standard finite difference method with explicit time-stepping
- **Boundary Conditions**: Supports both Dirichlet (fixed temperatures) and Neumann (zero-gradient) boundary conditions
- **Grid Update**: Updates entire grid at each time step sequentially
- **Output Location**: Saves results to `csv/seq/` directory
- **Best for**: Small to medium grid sizes (e.g., 50×50 to 256×256)

### Parallel Version ([parallele/parallel.cpp](parallele/parallel.cpp))

- **Method**: Domain decomposition using MPI
- **Parallelization Strategy**:
  - 2D Cartesian topology for process distribution
  - Each MPI process handles a subdomain of the grid
  - Ghost cell exchange between neighboring processes
- **Boundary Conditions**: Supports both Dirichlet (fixed temperatures) and Neumann (zero-gradient) boundary conditions
- **Communication**: Non-blocking sends/receives for ghost cell updates
- **Output Location**: Saves results to `csv/par/` directory
- **Best for**: Large grid sizes (e.g., 512×512 and above)

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
BOUNDARY_TYPE            # DIRICHLET or NEUMANN (both versions support both types)
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

### Boundary Condition Types

- **DIRICHLET**: Fixed temperature at boundaries (boundaries remain constant throughout simulation)
- **NEUMANN**: Zero-gradient at boundaries (no heat flux across boundaries, ∂T/∂n = 0)

## Usage

### 1. Generate Initial Conditions

**For custom large-scale simulations:**

```bash
python generate_ic_large.py
```

This creates a `initial_conditions_large.txt` file with a Gaussian temperature distribution:

- Default: 200×200 grid
- Hot spot at center with Gaussian decay
- Parameters: Lx=Ly=1.0, α=0.1, T_final=0.5

**For testing (generate multiple test cases):**

```bash
# Generate Dirichlet boundary condition test cases
python tests/tests_d.py

# Generate Neumann boundary condition test cases
python tests/tests_n.py
```

These scripts create test files in `initial_data/data_tests_d/` and `initial_data/data_tests_n/` directories with various grid sizes for weak/strong scaling analysis.

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
./sequential.exe sequential/initial_conditions.txt output_name
```

Output will be saved to `csv/seq/output_name.csv`

**Parallel (with 4 processes):**

```bash
mpirun -np 4 ./parallel.exe parallele/initial_conditions_large.txt output_name
```

Output will be saved to `csv/par/output_name_<BOUNDARY_TYPE>.csv` (e.g., `output_name_DIRICHLET.csv`)

**Notes:**

- Sequential output: `csv/seq/{output_prefix}.csv`
- Parallel output: `csv/par/{output_prefix}_{BC_TYPE}.csv`
- The boundary condition type is automatically appended to parallel output filenames

### 4. Visualize Results

Generate heatmap visualizations of the final temperature field:

```bash
python visualization.py csv/seq/output_name.csv plots
python visualization.py csv/par/output_name_DIRICHLET.csv plots
```

The visualization script:

- Reads CSV output files from `csv/seq/` or `csv/par/` directories
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

### Basic Workflow

Complete workflow for running and visualizing a simulation:

```bash
# Generate large initial conditions
python generate_ic_large.py

# Compile both versions
g++ sequential/sequential.cpp -o sequential.exe -O3
mpic++ parallele/parallel.cpp -o parallel.exe -O3

# Run sequential version
./sequential.exe sequential/initial_conditions.txt seq_result

# Run parallel version with 4 processes
mpirun -np 4 ./parallel.exe parallele/initial_conditions_large.txt par_result

# Visualize both results
python visualization.py csv/seq/seq_result.csv plots
python visualization.py csv/par/par_result_DIRICHLET.csv plots
```

### Testing Workflow with Different Boundary Conditions

```bash
# Generate test cases
python tests/tests_d.py  # Dirichlet tests
python tests/tests_n.py  # Neumann tests

# Compile programs
g++ sequential/sequential.cpp -o sequential.exe -O3
mpic++ parallele/parallel.cpp -o parallel.exe -O3

# Run with Dirichlet boundaries (256x256 grid)
./sequential.exe initial_data/data_tests_d/initial_weak_P1_256x256.txt seq_dirichlet_256
mpirun -np 4 ./parallel.exe initial_data/data_tests_d/initial_weak_P4_512x512.txt par_dirichlet_512

# Run with Neumann boundaries (256x256 grid)
./sequential.exe initial_data/data_tests_n/initial_weak_P1_256x256.txt seq_neumann_256
mpirun -np 4 ./parallel.exe initial_data/data_tests_n/initial_weak_P4_512x512.txt par_neumann_512

# Visualize results
python visualization.py csv/seq/seq_dirichlet_256.csv plots
python visualization.py csv/par/par_dirichlet_512_DIRICHLET.csv plots
python visualization.py csv/seq/seq_neumann_256.csv plots
python visualization.py csv/par/par_neumann_512_NEUMANN.csv plots
```

## Output Files

### CSV Output Structure

- **Sequential**: Files saved to `csv/seq/{output_prefix}.csv`
- **Parallel**: Files saved to `csv/par/{output_prefix}_{BC_TYPE}.csv`
- **Format**: Comma-separated temperature values (Ny rows × Nx columns)
- **Data Layout**: Each row represents a horizontal line of the grid (j-index), each column is a vertical position (i-index)

### Visualization Output

- **Plot Images**: PNG files with timestamp saved to `plots/` directory
- **Naming Convention**: Based on input CSV filename with timestamp

## Features Summary

### Boundary Conditions

Both sequential and parallel implementations now support:

1. **Dirichlet Boundaries**: Fixed temperature at domain edges

   - Boundaries remain constant throughout simulation
   - Useful for modeling heat sources/sinks at boundaries

2. **Neumann Boundaries**: Zero-gradient condition (∂T/∂n = 0)
   - No heat flux across boundaries
   - Useful for insulated/isolated systems

### Automatic Output Organization

- Sequential results → `csv/seq/`
- Parallel results → `csv/par/`
- Visualizations → `plots/`

### Test Data Generation

- `tests/tests_d.py`: Generates Dirichlet test cases in `initial_data/data_tests_d/`
- `tests/tests_n.py`: Generates Neumann test cases in `initial_data/data_tests_n/`
- Includes reference cases and weak/strong scaling test configurations

## Notes

- Both sequential and parallel versions support Dirichlet and Neumann boundary conditions
- Initial conditions can be customized by modifying `generate_ic_large.py` or test scripts
- For best parallel performance, use grid sizes that are evenly divisible by the number of processes
- The visualization colormap can be changed in `visualization.py` (currently uses "viridis")
- Automatic time-step adjustment ensures numerical stability based on CFL condition
