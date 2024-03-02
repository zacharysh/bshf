#include <iostream>

#include <vector>

#include <cmath> // sqrt

#include <fstream>

#include "IO/io.hpp"

#include "wavefunction/wavefunction.hpp"
#include "wavefunction/atom.hpp"
#include "wavefunction/solve.hpp"

#include "math/matrix.hpp"
#include "math/math_other.hpp"

int main(int argc, char **argv)
{
    /*
    if (argc != 4)
    {
        std::cerr << "Incorrect number of args. Requires 3: N, tmax, v0\n";
        return EXIT_FAILURE;
    }
    */
    
    const double r_min = 1.0e-8;
    const double r_max = 20.0;
    const int k_spline = 7;
    const int n_spline = 63;
    const int N_grid = 10000;

    std::cout << "> Constructing linear grid: r0 = " << r_min << ", r_max = " << r_max << ", N_grid = " << N_grid << "... ";
    std::vector<double> r_grid = construct_grid_linear(r_min, r_max, N_grid);
    std::cout << IO::done;

    SplineBasis basis(r_grid, k_spline, n_spline);
    Atom::AtomicSystem Li(3, r_grid);
    Wavefunction psi(1, 0, 0); // H-like 1s orbital

    solve_hydrogen_like(Li, psi, basis);

    
    std::ofstream ofs;
    ofs.open("output/Li_1s.txt", std::ofstream::out | std::ofstream::app);

    //ofs << "r, psi(r)\n";

    
    for(int i = 0; i < N_grid; ++i)
    {
        ofs << r_grid.at(i) << ", " << psi.P.at(i) << ", " << psi.amplitude.at(i) << '\n';
    }

    ofs.close();

    
    return EXIT_SUCCESS;
}