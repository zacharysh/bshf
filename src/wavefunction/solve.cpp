#include "solve.hpp"


auto constructHamiltonian(Atom::AtomicSystem system, Wavefunction psi, SplineBasis *basis) -> SquareMatrix<double> *
{
    std::cout << "Constructing Hamiltonian... ";
    int N_b = basis->num_spl;

    // fix me up
    // generate centrifugal term
    std::vector<double> centrifugal{};//(basis->r_grid.size());
    centrifugal.reserve(basis->r_grid.size());
    for (auto r : basis->r_grid)
        centrifugal.push_back(psi.l * (psi.l + 1)/ (2 * r * r));

    SquareMatrix<double> *H = new SquareMatrix<double>(N_b);

    for(int i = 0; i < N_b; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            (*H)(i, j) = 
                  integrate_trap(basis->r_grid, basis->bspl_derivative.at(i) * basis->bspl_derivative.at(j)) / 2.0;
                + integrate_trap(basis->r_grid, basis->bspl.at(i) * system.potential * basis->bspl.at(j))
                + integrate_trap(basis->r_grid, basis->bspl.at(i) * centrifugal * basis->bspl.at(j));
        }
    }
    std::cout << "done.\n";
    return H;
}

auto constructBMatrix(SplineBasis *basis) -> SquareMatrix<double>*
{
    std::cout << "Constructing BMatrix... ";
    auto N_b = basis->num_spl;

    SquareMatrix<double> *B = new SquareMatrix<double>(N_b);

    for(int i = 0; i < N_b; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            (*B)(i,j) = integrate_trap(basis->r_grid, basis->bspl.at(i) * basis->bspl.at(j));
        }
    }

    std::cout << "done.\n";
    return B;
}

auto solveHydrogenlike(Wavefunction *psi, SquareMatrix<double> *Hamiltonian, SquareMatrix<double> *BMatrix) -> void
{
    double *energies = MatrixTools::solveEigenSystem(Hamiltonian, BMatrix);

    int Z = 3;

    size_t N = Hamiltonian->get_size();

    std::cout << "--- Energies ---" << '\n';
    std::cout << " n | En    | predicted\n"; 
    for (int i=1; i <= 10; i++)
    {
        printf("%2d | %1.3f | %1.3f", i, energies[i-1], Z * Z / (2.0 * i * i));
        std::cout << "\n";
    }

    delete [] energies;
}