#include <iostream>

#include <vector>

#include <fstream>

#include "wavefunction/wavefunction.hpp"
#include "wavefunction/atom.hpp"
#include "wavefunction/solve.hpp"

#include "math/matrix.hpp"

int main(int argc, char **argv)
{
    /*
    if (argc != 4)
    {
        std::cerr << "Incorrect number of args. Requires 3: N, tmax, v0\n";
        return EXIT_FAILURE;
    }
    */
    
    
    double r_min = 1.0e-5;
    double r_max = 20.0;
    int k_spline = 7;
    int n_spline = 60;

    int N_grid = 2000;

    std::cout << "Constructing linear grid: r0 = " << r_min << ", r_max = " << r_max << ", N_grid = " << N_grid << "... ";
    std::vector<double> r_grid = constructGridLin(r_min, r_max, N_grid);
    std::cout << "done.\n";

    std::cout << "Constructing BSpline basis, k = " << k_spline << ", n = " << n_spline << "... ";
    SplineBasis basis(r_grid, k_spline, n_spline);
    std::cout << "done.\n";
    
    Atom::AtomicSystem Li(3, r_grid);

    Wavefunction psi(1, 0, 0); // H-like 1s orbital

    SquareMatrix<double> *Hamiltonian = constructHamiltonian(Li, psi, &basis);
    SquareMatrix<double> *BMatrix = constructBMatrix(&basis);

    solveHydrogenlike(&psi, Hamiltonian, BMatrix);

    delete Hamiltonian;
    delete BMatrix;

    /*
    std::cout << "--- Energies ---" << '\n';
    for (size_t i=0; i < psi.energy.size(); i++)
    {
        std::cout << "( " << i << ", " << psi.energy.at(i) << " )\n";
    }
    std::cout << '\n';
    */

    return EXIT_SUCCESS;
}