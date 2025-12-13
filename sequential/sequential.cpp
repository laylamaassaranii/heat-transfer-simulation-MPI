#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <limits>   // for numeric_limits

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

int main(int argc, char **argv)
{
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " initial_conditions.txt output_prefix\n";
        return 1;
    }

    cout << "Entered the function main\n";

    // --- Simple switch to enable/disable debug output ---
    const bool DEBUG = true;
    const int DEBUG_PRINT_EVERY = 1000; // print every 1000 time steps

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

    if (DEBUG) {
        cout << "[DEBUG] Header read from file:\n";
        cout << "        Nx = " << Nx << ", Ny = " << Ny << "\n";
        cout << "        Lx = " << Lx << ", Ly = " << Ly << "\n";
        cout << "        alpha = " << alpha << "\n";
        cout << "        T_final = " << T_final << "\n";
        cout << "        dt_in = " << dt_in << "\n";
        cout << "        bc_str = " << bc_str << "\n";
    }

    if (Nx < 3 || Ny < 3) {
        cerr << "Error: Nx and Ny must be >= 3 (need interior points).\n";
        return 1;
    }

    Grid T(Nx, vector<double>(Ny));
    Grid T_new(Nx, vector<double>(Ny));

    int count_vals = 0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (!(in >> T[i][j])) {
                cerr << "Error: not enough temperature values in file.\n";
                cerr << "       Only read " << count_vals << " values out of Nx*Ny = " 
                     << Nx * Ny << "\n";
                return 1;
            }
            ++count_vals;
        }
    }

    if (DEBUG) {
        cout << "[DEBUG] Successfully read " << count_vals 
             << " temperature values (expected " << Nx * Ny << ")\n";
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

    // Save Dirichlet boundary values
    vector<double> left(Ny), right(Ny), bottom(Nx), top(Nx);

    for (int j = 0; j < Ny; ++j) {
        left[j]  = T[0][j];
        right[j] = T[Nx - 1][j];
    }
    for (int i = 0; i < Nx; ++i) {
        bottom[i] = T[i][0];
        top[i]    = T[i][Ny - 1];
    }

    if (DEBUG) {
        cout << "[DEBUG] Example initial values:\n";
        cout << "        T(center)   = T[" << Nx/2 << "][" << Ny/2 << "] = "
             << T[Nx/2][Ny/2] << "\n";
        cout << "        T(left mid) = T[0][" << Ny/2 << "] = " 
             << T[0][Ny/2] << " (Dirichlet)\n";
        cout << "        T(right mid)= T[" << Nx-1 << "][" << Ny/2 << "] = " 
             << T[Nx-1][Ny/2] << " (Dirichlet)\n";
    }

    // --- Time stepping ---
    for (int n = 0; n < Nt; ++n) {
        apply_dirichlet(T, left, right, bottom, top);

        // Optional: compute min/max to ensure solution stays finite
        double minT = numeric_limits<double>::infinity();
        double maxT = -numeric_limits<double>::infinity();

        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double Tij   = T[i][j];
                double Tip1j = T[i + 1][j];
                double Tim1j = T[i - 1][j];
                double Tijp1 = T[i][j + 1];
                double Tijm1 = T[i][j - 1];

                double lap =
                    (Tip1j - 2.0 * Tij + Tim1j) / (dx * dx) +
                    (Tijp1 - 2.0 * Tij + Tijm1) / (dy * dy);

                double val = Tij + alpha * dt * lap;
                T_new[i][j] = val;

                if (val < minT) minT = val;
                if (val > maxT) maxT = val;
            }
        }

        apply_dirichlet(T_new, left, right, bottom, top);

        T.swap(T_new);

        // Debug print every DEBUG_PRINT_EVERY iterations
        if (DEBUG && (n % DEBUG_PRINT_EVERY == 0)) {
            cout << "[DEBUG] Step n = " << n
                 << ", time t = " << n*dt
                 << ", min(T) = " << minT
                 << ", max(T) = " << maxT << "\n";

            // Check that Dirichlet boundaries are still enforced
            cout << "        T(0,Ny/2)   = " << T[0][Ny/2]
                 << " (should equal left[Ny/2] = " << left[Ny/2] << ")\n";
            cout << "        T(Nx-1,Ny/2)= " << T[Nx-1][Ny/2]
                 << " (should equal right[Ny/2] = " << right[Ny/2] << ")\n";
        }
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
