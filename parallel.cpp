#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

int local_size(int N, int coord, int dim) {
    int base = N / dim;
    int rem  = N % dim;
    return base + (coord < rem ? 1 : 0);
}

int global_start(int N, int coord, int dim) {
    int base = N / dim;
    int rem  = N % dim;
    if (coord < rem) {
        return coord * (base + 1);
    } else {
        return rem * (base + 1) + (coord - rem) * base;
    }
}

inline int idx(int i, int j, int local_ny_with_halo) {
    return i * local_ny_with_halo + j;
}

void apply_neumann_ghosts(vector<double> &T, int local_nx, int local_ny, int Ny_with_halo, int coords[2], int dims[2]) {
    if (coords[0] == 0) {
        for (int j = 1; j <= local_ny; ++j) {
            T[idx(0, j, Ny_with_halo)] = T[idx(1, j, Ny_with_halo)];
        }
    }
    if (coords[0] == dims[0] - 1) {
        for (int j = 1; j <= local_ny; ++j) {
            T[idx(local_nx + 1, j, Ny_with_halo)] = T[idx(local_nx, j, Ny_with_halo)];
        }
    }
    if (coords[1] == 0) {
        for (int i = 1; i <= local_nx; ++i) {
            T[idx(i, 0, Ny_with_halo)] = T[idx(i, 1, Ny_with_halo)];
        }
    }
    if (coords[1] == dims[1] - 1) {
        for (int i = 1; i <= local_nx; ++i) {
            T[idx(i, local_ny + 1, Ny_with_halo)] = T[idx(i, local_ny, Ny_with_halo)];
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm world = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(world, &rank);
    MPI_Comm_size(world, &size);

    if (rank == 0 && argc < 3) {
        cerr << "Usage: mpirun -np <P> " << argv[0] << " initial_conditions_large.txt output_prefix\n";
        MPI_Abort(world, 1);
    }

    string input_path;
    string output_prefix;
    if (rank == 0) {
        input_path    = argv[1];
        output_prefix = argv[2];
    }

    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(world, 2, dims, periods, 1, &cart_comm);

    int cart_rank;
    int coords[2];
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

    int nbr_left, nbr_right, nbr_down, nbr_up;
    MPI_Cart_shift(cart_comm, 0, 1, &nbr_left,  &nbr_right);
    MPI_Cart_shift(cart_comm, 1, 1, &nbr_down, &nbr_up);

    int Nx = 0, Ny = 0;
    double Lx = 0.0, Ly = 0.0;
    double alpha = 0.0;
    double T_final = 0.0;
    double dt_in = 0.0;
    int bc_type = 0;

    vector<double> global_T;

    if (cart_rank == 0) {
        ifstream in(input_path);
        if (!in) {
            cerr << "Error: cannot open input file " << input_path << "\n";
            MPI_Abort(cart_comm, 1);
        }

        string bc_str;
        in >> Nx >> Ny;
        in >> Lx >> Ly;
        in >> alpha;
        in >> T_final;
        in >> dt_in;
        in >> bc_str;

        if (!in) {
            cerr << "Error: invalid header in input file.\n";
            MPI_Abort(cart_comm, 1);
        }

        transform(bc_str.begin(), bc_str.end(), bc_str.begin(), ::toupper);
        if (bc_str == "NEUMANN") {
            bc_type = 1;
        } else {
            bc_type = 0;
        }

        global_T.resize(Nx * Ny);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (!(in >> global_T[j * Nx + i])) {
                    cerr << "Error: not enough temperature values in file.\n";
                    MPI_Abort(cart_comm, 1);
                }
            }
        }
        in.close();
    }

    MPI_Bcast(&Nx, 1, MPI_INT, 0, cart_comm);
    MPI_Bcast(&Ny, 1, MPI_INT,0 , cart_comm);
    MPI_Bcast(&Lx, 1 , MPI_DOUBLE, 0, cart_comm);
    MPI_Bcast(&Ly, 1, MPI_DOUBLE, 0, cart_comm);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, cart_comm);
    MPI_Bcast(&T_final, 1, MPI_DOUBLE, 0, cart_comm);
    MPI_Bcast(&dt_in, 1, MPI_DOUBLE, 0, cart_comm);
    MPI_Bcast(&bc_type, 1, MPI_INT, 0, cart_comm);

    if (Nx < 3 || Ny < 3) {
        if (cart_rank == 0) {
            cerr << "Error: Nx and Ny must be >= 3.\n";
        }
        MPI_Abort(cart_comm, 1);
    }

    int local_nx = local_size(Nx, coords[0], dims[0]);
    int local_ny = local_size(Ny, coords[1], dims[1]);

    int global_i_start = global_start(Nx, coords[0], dims[0]);
    int global_j_start = global_start(Ny, coords[1], dims[1]);

    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);

    double denom = (1.0 / (dx * dx) + 1.0 / (dy * dy));
    double dt_max = 0.5 / (alpha * denom);

    double dt;
    if (dt_in <= 0.0 || dt_in > dt_max) {
        dt = 0.9 * dt_max;
    } else {
        dt = dt_in;
    }
    int Nt = static_cast<int>(T_final / dt);
    double T_final_effective = Nt * dt;

    if (cart_rank == 0) {
        cout << "[INFO] Global Nx=" << Nx << ", Ny=" << Ny << ", Lx=" << Lx << ", Ly=" << Ly << "\n";
        cout << "[INFO] alpha=" << alpha << ", dt=" << dt << " (dt_max=" << dt_max << "), Nt=" << Nt << ", T_final_effective=" << T_final_effective << "\n";
        cout << "[INFO] dims = [" << dims[0] << " x " << dims[1] << "]\n";
    }

    if (Nt <= 0) {
        if (cart_rank == 0) {
            cerr << "Error: Nt <= 0 (T_final too small or dt too large).\n";
        }
        MPI_Abort(cart_comm, 1);
    }

    if (cart_rank == 0) {
        cout << "[INFO] Rank 0 will distribute the initial field...\n";
    }

    int Ny_with_halo = local_ny + 2;
    int Nx_with_halo = local_nx + 2;
    vector<double> T   (Nx_with_halo * Ny_with_halo, 0.0);
    vector<double> Tnew(Nx_with_halo * Ny_with_halo, 0.0);

    vector<double> local_block(local_nx * local_ny);

    if (cart_rank == 0) {
        for (int p = 0; p < size; ++p) {
            int c[2];
            MPI_Cart_coords(cart_comm, p, 2, c);

            int ln_x = local_size(Nx, c[0], dims[0]);
            int ln_y = local_size(Ny, c[1], dims[1]);
            int gi0  = global_start(Nx, c[0], dims[0]);
            int gj0  = global_start(Ny, c[1], dims[1]);

            vector<double> sendbuf(ln_x * ln_y);
            for (int lj = 0; lj < ln_y; ++lj) {
                for (int li = 0; li < ln_x; ++li) {
                    int gi = gi0 + li;
                    int gj = gj0 + lj;
                    sendbuf[lj * ln_x + li] = global_T[gj * Nx + gi];
                }
            }

            if (p == 0) {
                local_block = sendbuf;
            } else {
                MPI_Send(sendbuf.data(), ln_x * ln_y, MPI_DOUBLE, p, 100, cart_comm);
            }
        }
    } else {
        MPI_Recv(local_block.data(), local_nx * local_ny, MPI_DOUBLE,
                 0, 100, cart_comm, MPI_STATUS_IGNORE);
    }

    for (int lj = 0; lj < local_ny; ++lj) {
        for (int li = 0; li < local_nx; ++li) {
            int I = li + 1;
            int J = lj + 1;
            T[idx(I, J, Ny_with_halo)] = local_block[lj * local_nx + li];
        }
    }

    vector<double> send_left (local_ny), recv_left (local_ny);
    vector<double> send_right(local_ny), recv_right(local_ny);
    vector<double> send_down(local_nx), recv_down(local_nx);
    vector<double> send_up  (local_nx), recv_up  (local_nx);

    for (int n = 0; n < Nt; ++n) {
        MPI_Request reqs[8];
        int req_count = 0;

        if (nbr_left != MPI_PROC_NULL) {
            for (int j = 1; j <= local_ny; ++j) {
                send_left[j - 1] = T[idx(1, j, Ny_with_halo)];
            }
            MPI_Irecv(recv_left.data(), local_ny, MPI_DOUBLE,
                      nbr_left,  10, cart_comm, &reqs[req_count++]);
            MPI_Isend(send_left.data(), local_ny, MPI_DOUBLE,
                      nbr_left,  11, cart_comm, &reqs[req_count++]);
        }
        if (nbr_right != MPI_PROC_NULL) {
            for (int j = 1; j <= local_ny; ++j) {
                send_right[j - 1] = T[idx(local_nx, j, Ny_with_halo)];
            }
            MPI_Irecv(recv_right.data(), local_ny, MPI_DOUBLE,
                      nbr_right, 11, cart_comm, &reqs[req_count++]);
            MPI_Isend(send_right.data(), local_ny, MPI_DOUBLE,
                      nbr_right, 10, cart_comm, &reqs[req_count++]);
        }

        if (nbr_down != MPI_PROC_NULL) {
            for (int i = 1; i <= local_nx; ++i) {
                send_down[i - 1] = T[idx(i, 1, Ny_with_halo)];
            }
            MPI_Irecv(recv_down.data(), local_nx, MPI_DOUBLE,
                      nbr_down, 20, cart_comm, &reqs[req_count++]);
            MPI_Isend(send_down.data(), local_nx, MPI_DOUBLE,
                      nbr_down, 21, cart_comm, &reqs[req_count++]);
        }
        if (nbr_up != MPI_PROC_NULL) {
            for (int i = 1; i <= local_nx; ++i) {
                send_up[i - 1] = T[idx(i, local_ny, Ny_with_halo)];
            }
            MPI_Irecv(recv_up.data(), local_nx, MPI_DOUBLE,
                      nbr_up, 21, cart_comm, &reqs[req_count++]);
            MPI_Isend(send_up.data(), local_nx, MPI_DOUBLE,
                      nbr_up, 20, cart_comm, &reqs[req_count++]);
        }

        if (req_count > 0) {
            MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);
        }

        if (nbr_left != MPI_PROC_NULL) {
            for (int j = 1; j <= local_ny; ++j) {
                T[idx(0, j, Ny_with_halo)] = recv_left[j - 1];
            }
        }
        if (nbr_right != MPI_PROC_NULL) {
            for (int j = 1; j <= local_ny; ++j) {
                T[idx(local_nx + 1, j, Ny_with_halo)] = recv_right[j - 1];
            }
        }
        if (nbr_down != MPI_PROC_NULL) {
            for (int i = 1; i <= local_nx; ++i) {
                T[idx(i, 0, Ny_with_halo)] = recv_down[i - 1];
            }
        }
        if (nbr_up != MPI_PROC_NULL) {
            for (int i = 1; i <= local_nx; ++i) {
                T[idx(i, local_ny + 1, Ny_with_halo)] = recv_up[i - 1];
            }
        }

        if (bc_type == 1) {
            apply_neumann_ghosts(T, local_nx, local_ny, Ny_with_halo, coords, dims);
        }

        for (int i = 1; i <= local_nx; ++i) {
            int gi = global_i_start + (i - 1);
            for (int j = 1; j <= local_ny; ++j) {
                int gj = global_j_start + (j - 1);

                double Tij = T[idx(i, j, Ny_with_halo)];

                if (bc_type == 0 &&
                    (gi == 0 || gi == Nx - 1 || gj == 0 || gj == Ny - 1)) {
                    Tnew[idx(i, j, Ny_with_halo)] = Tij;
                    continue;
                }

                double Tip1j = T[idx(i + 1, j, Ny_with_halo)];
                double Tim1j = T[idx(i - 1, j, Ny_with_halo)];
                double Tijp1 = T[idx(i, j + 1, Ny_with_halo)];
                double Tijm1 = T[idx(i, j - 1, Ny_with_halo)];

                double lap =
                    (Tip1j - 2.0 * Tij + Tim1j) / (dx * dx) +
                    (Tijp1 - 2.0 * Tij + Tijm1) / (dy * dy);

                Tnew[idx(i, j, Ny_with_halo)] = Tij + alpha * dt * lap;
            }
        }

        swap(T, Tnew);
    }

    for (int lj = 0; lj < local_ny; ++lj) {
        for (int li = 0; li < local_nx; ++li) {
            int I = li + 1;
            int J = lj + 1;
            local_block[lj * local_nx + li] = T[idx(I, J, Ny_with_halo)];
        }
    }

    if (cart_rank == 0) {
        vector<double> global_final(Nx * Ny);

        for (int lj = 0; lj < local_ny; ++lj) {
            for (int li = 0; li < local_nx; ++li) {
                int gi = global_i_start + li;
                int gj = global_j_start + lj;
                global_final[gj * Nx + gi] = local_block[lj * local_nx + li];
            }
        }

        for (int p = 1; p < size; ++p) {
            int c[2];
            MPI_Cart_coords(cart_comm, p, 2, c);
            int ln_x = local_size(Nx, c[0], dims[0]);
            int ln_y = local_size(Ny, c[1], dims[1]);
            int gi0  = global_start(Nx, c[0], dims[0]);
            int gj0  = global_start(Ny, c[1], dims[1]);

            vector<double> recv_block(ln_x * ln_y);
            MPI_Recv(recv_block.data(), ln_x * ln_y, MPI_DOUBLE,
                     p, 200, cart_comm, MPI_STATUS_IGNORE);

            for (int lj = 0; lj < ln_y; ++lj) {
                for (int li = 0; li < ln_x; ++li) {
                    int gi = gi0 + li;
                    int gj = gj0 + lj;
                    global_final[gj * Nx + gi] = recv_block[lj * ln_x + li];
                }
            }
        }

        string out_path = output_prefix + "_final.csv";
        ofstream out(out_path);
        if (!out) {
            cerr << "Error: cannot open output file " << out_path << "\n";
        } else {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    out << global_final[j * Nx + i];
                    if (i < Nx - 1) out << ',';
                }
                out << '\n';
            }
            out.close();
            cout << "[INFO] Final field written to " << out_path << "\n";
        }
    } else {
        MPI_Send(local_block.data(), local_nx * local_ny, MPI_DOUBLE,
                 0, 200, cart_comm);
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
