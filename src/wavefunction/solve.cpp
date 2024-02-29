#include "solve.hpp"


auto constructHamiltonian(Atom::AtomicSystem system, Wavefunction psi, SplineBasis &basis) -> SquareMatrix<double>
{
    std::cout << "  > Constructing Hamiltonian... ";
    int N_b = basis.num_spl;

    // fix me up
    // generate centrifugal term
    std::vector<double> centrifugal{};//(basis->r_grid.size());
    centrifugal.reserve(basis.r_grid.size());
    for (auto r : basis.r_grid)
        centrifugal.push_back(psi.l * (psi.l + 1)/ (2 * r * r));

    SquareMatrix<double> H(N_b);

    #pragma omp parallel for
    for(int i = 0; i < N_b; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) = 
                  integrate_trap(basis.r_grid, basis.bspl_derivative.at(i) * basis.bspl_derivative.at(j)) / 2.0;
                + integrate_trap(basis.r_grid, basis.bspl.at(i) * system.potential * basis.bspl.at(j))
                + integrate_trap(basis.r_grid, basis.bspl.at(i) * centrifugal * basis.bspl.at(j));
        }
    }
    std::cout << "done.\n";
    return H;
}

auto constructBMatrix(SplineBasis &basis) -> SquareMatrix<double>
{
    std::cout << "  > Constructing BMatrix... ";
    auto N_b = basis.num_spl;

    SquareMatrix<double> B(N_b);

    #pragma omp parallel for
    for(int i = 0; i < N_b; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            B(i,j) = integrate_trap(basis.r_grid, basis.bspl.at(i) * basis.bspl.at(j));
        }
    }

    std::cout << "done.\n";
    return B;
}

auto solveHydrogenlike(Atom::AtomicSystem &atom, Wavefunction &psi, SplineBasis &basis) -> void
{
    auto Hamiltonian = constructHamiltonian(atom, psi, basis);
    auto BMatrix = constructBMatrix(basis);

    auto energies = MatrixTools::solveEigenSystem(Hamiltonian, BMatrix);

    std::size_t N = Hamiltonian.get_size();

    std::cout << "--- Energies ---" << '\n';
    std::cout << " n | En    | predicted\n"; 
    for (int i=1; i <= 10; ++i)
    {
        printf("%2d | %1.3f | %1.3f", i, energies.at(N-i), atom.Z * atom.Z / (2.0 * i * i));
        std::cout << "\n";
    }

    for (std::size_t i = 0; i < N; ++i)
        psi.expansion_coeffs.push_back(Hamiltonian(psi.n - 1, N - i));

}