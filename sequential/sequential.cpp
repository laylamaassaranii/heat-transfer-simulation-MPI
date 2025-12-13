#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <limits>

using namespace std;
using Grid = vector<vector<double>>;

void apply_dirichlet(Grid &T,
                     const vector<double> &left,
                     const vector<double> &right,
                     const vector<double> &bottom,
                     const vector<double> &top)
{
    int Nx = static_cast<int>(T.size());
    int Ny = static_cast<int>(T[0].size());

    for (int j = 0; j < Ny; ++j) {
        T[0][j] = left[j];
        T[Nx - 1][j] = right[j];
    }

    for (int i = 0; i < Nx; ++i) {
        T[i][0] = bottom[i];
        T[i][Ny - 1] = top[i];
    }
}

void apply_neumann(Grid &T)
{
    int Nx = static_cast<int>(T.size());
    int Ny = static_cast<int>(T[0].size());

    for (int j = 0; j < Ny; ++j) {
        T[0][j]      = T[1][j];
        T[Nx - 1][j] = T[Nx - 2][j];
    }

    for (int i = 0; i < Nx; ++i) {
        T[i][0]      = T[i][1];
        T[i][Ny - 1] = T[i][Ny - 2];
    }
}

int main(int argc, char **argv)
    {
        if (argc < 3) {
            cerr << "Usage: " << argv[0] << " initial_conditions.txt output_prefix\n";
            return 1;
        }

        string input_path = argv[1];
        string output_pref = argv[2];

        ifstream in(input_path);
        if (!in) {
            cerr << "Error: cannot open input file " << input_path << "\n";
            return 1;
        }

        int Nx, Ny;
        double Lx, Ly;
        double alpha;
        double T_final;
        double dt_in;
        string bc_str;

        in >> Nx >> Ny;
        in >> Lx >> Ly;
        in >> alpha;
        in >> T_final;
        in >> dt_in;
        in >> bc_str;

        if (!in) {
            cerr << "Error: invalid header in input file.\n";
            return 1;
        }

        if (Nx < 3 || Ny < 3) {
            cerr << "Error: Nx and Ny must be >= 3 (need interior points).\n";
            return 1;
        }

        int bc_type = -1;
        string bc_upper = bc_str;
        transform(bc_upper.begin(), bc_upper.end(), bc_upper.begin(), ::toupper);

        if (bc_upper == "DIRICHLET") {
            bc_type = 0;
            cout << "Using DIRICHLET boundary conditions (bc_type = 0).\n";
        }
        else if (bc_upper == "NEUMANN") {
            bc_type = 1;
            cout << "Using NEUMANN boundary conditions (bc_type = 1).\n";
        }
        else {
            cerr << "Error: boundary condition '" << bc_str
                << "' is not implemented. Supported BCs are: DIRICHLET, NEUMANN.\n";
            return 1;
        }

        Grid T(Nx, vector<double>(Ny));
        Grid T_new(Nx, vector<double>(Ny));

        int count_vals = 0;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (!(in >> T[i][j])) {
                    cerr << "Error: not enough temperature values in file.\n";
                    return 1;
                }
                ++count_vals;
            }
        }

        in.close();

        double dx = Lx / (Nx - 1);
        double dy = Ly / (Ny - 1);

        double denom = (1.0 / (dx * dx) + 1.0 / (dy * dy));
        double dt_max = 0.5 / (alpha * denom);

        double dt;
        if (dt_in <= 0.0 || dt_in > dt_max) {
            dt = 0.9 * dt_max;
            cout << "Adjusting dt to stable value: dt = " << dt << " (dt_max = " << dt_max << ")\n";
        } else {
            dt = dt_in;
        }

        int Nt = static_cast<int>(T_final / dt);
        double T_final_effective = Nt * dt;

        cout << "Nx=" << Nx << ", Ny=" << Ny << ", dx=" << dx << ", dy=" << dy << "\n";
        cout << "alpha=" << alpha << ", dt=" << dt << ", Nt=" << Nt
            << ", T_final_effective=" << T_final_effective << "\n";

        if (Nt <= 0) {
            cerr << "Error: Nt <= 0 (T_final too small or dt too large).\n";
            return 1;
        }

        cout << "Nx=" << Nx << ", Ny=" << Ny
            << ", dx=" << dx << ", dy=" << dy << "\n";
        cout << "alpha=" << alpha << ", dt=" << dt
            << ", Nt=" << Nt
            << ", T_final_effective=" << T_final_effective << "\n";

        vector<double> left, right, bottom, top;
        if (bc_type == 0) {
            left.assign(Ny, 0.0);
            right.assign(Ny, 0.0);
            bottom.assign(Nx, 0.0);
            top.assign(Nx, 0.0);

        for (int j = 0; j < Ny; ++j) {
            left[j]  = T[0][j];
            right[j] = T[Nx - 1][j];
        }
        for (int i = 0; i < Nx; ++i) {
            bottom[i] = T[i][0];
            top[i]    = T[i][Ny - 1];
        }

        for (int n = 0; n < Nt; ++n) {
            if (bc_type == 0) {
                apply_dirichlet(T, left, right, bottom, top);
            } else if (bc_type == 1) {
                apply_neumann(T);
            }

            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    double Tij = T[i][j];
                    double Tip1j = T[i + 1][j];
                    double Tim1j = T[i - 1][j];
                    double Tijp1 = T[i][j + 1];
                    double Tijm1 = T[i][j - 1];

                    double lap =
                        (Tip1j - 2.0 * Tij + Tim1j) / (dx * dx) +
                        (Tijp1 - 2.0 * Tij + Tijm1) / (dy * dy);

                    T_new[i][j] = Tij + alpha * dt * lap;
                }
            }

            if (bc_type == 0) {
                apply_dirichlet(T_new, left, right, bottom, top);
            } else {
                apply_neumann(T_new);
            }

            T.swap(T_new);
        }

        string out_path = output_pref + ".csv";
        ofstream out(out_path);
        if (!out) {
            cerr << "Error: cannot open output file " << out_path << "\n";
            return 1;
        }

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                out << T[i][j];
                if (i < Nx - 1)
                    out << ',';
            }
            out << '\n';
        }

        out.close();
        cout << "Final field written to " << out_path << "\n";
        return 0;
    }
}
