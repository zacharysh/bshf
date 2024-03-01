#include "solve.hpp"


auto constructHamiltonian(Atom::AtomicSystem system, Wavefunction psi, SplineBasis &basis) -> SquareMatrix<double>
{
    std::cout << "  > Constructing Hamiltonian... ";

    // fix me up
    // generate centrifugal term
    std::vector<double> centrifugal(basis.r_grid.size());

    for (std::size_t i = 0; i < basis.r_grid.size(); ++i) 
        centrifugal.at(i) =  psi.l * (psi.l + 1)  / (2.0 * basis.r_grid.at(i) * basis.r_grid.at(i));
    

    SquareMatrix<double> H(basis.num_spl);

    #pragma omp parallel for
    for(int i = 0; i < basis.num_spl; ++i)
    {
        for(int j = 0; j < basis.num_spl; ++j) // <= i
        {
            H(i, j) = simpson_linear(basis.dr, basis.bspl_derivative.at(i) * basis.bspl_derivative.at(j)) / 2.0;
                    + simpson_linear(basis.dr, basis.bspl.at(i) * system.potential * basis.bspl.at(j))
                    + simpson_linear(basis.dr, basis.bspl.at(i) * centrifugal * basis.bspl.at(j));
        }
    }
    std::cout << "done.\n";
    return H;
}

auto constructBMatrix(SplineBasis &basis) -> SquareMatrix<double>
{
    std::cout << "  > Constructing BMatrix... ";

    SquareMatrix<double> B(basis.num_spl);

    #pragma omp parallel for
    for(int i = 0; i < basis.num_spl; ++i)
    {
        for(int j = 0; j < basis.num_spl; ++j)
        {
            B(i,j) = simpson_linear(basis.dr, basis.bspl.at(i) * basis.bspl.at(j));
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

    std::cout << "\n\n";

    std::cout << "--- Energies ---" << '\n';
    std::cout << " n | En    | predicted\n"; 
    for (int i=1; i <= 10; ++i)
    {
        printf("%2d | %1.3f | %1.3f", i, energies.at(i-1), atom.Z * atom.Z / (2.0 * i * i));
        std::cout << "\n";
    }

    for (std::size_t i = 0; i < N; ++i)
    {
        psi.expansion_coeffs.push_back(Hamiltonian(psi.n - 1, i));
    }

    psi.amplitude.resize(basis.r_grid.size());


    //#pragma omp parallel for
    for(int i = 0; i < basis.num_spl; ++i)
    {
        for(std::size_t j = 0; j < basis.r_grid.size(); ++j)
        {
            psi.amplitude.at(j) += basis.bspl.at(i).at(j) * psi.expansion_coeffs.at(i);
        }        
    }

    std::cout << "\n\u222Bdr|P(r)|\u00B2 = " << simpson_linear(basis.dr, psi.amplitude * psi.amplitude) << ".\n";

}