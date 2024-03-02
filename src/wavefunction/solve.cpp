#include "solve.hpp"


auto construct_hamiltonian_matrix(Atom::AtomicSystem system, Wavefunction psi, const SplineBasis &basis) -> SquareMatrix<double>
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
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) = trapz_linear(basis.dr, basis.bspl_derivative.at(i) * basis.bspl_derivative.at(j)) / 2.0 // IT WAS A FUCKING SEMICOLON
                    + trapz_linear(basis.dr, basis.bspl.at(i) * system.potential * basis.bspl.at(j))
                    + trapz_linear(basis.dr, basis.bspl.at(i) * centrifugal * basis.bspl.at(j));

            // Only need to calculate bottom half.
            H(j,i) = H(i,j);
        }
    }
    std::cout << IO::done;
    return H;
}

auto construct_spline_matrix(const SplineBasis &basis) -> SquareMatrix<double>
{
    std::cout << "  > Constructing BMatrix... ";

    SquareMatrix<double> B(basis.num_spl);

    #pragma omp parallel for
    for(int i = 0; i < basis.num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            B(i,j) = trapz_linear(basis.dr, basis.bspl.at(i) * basis.bspl.at(j));

            // Only need to calculate bottom half.
            B(j,i) = B(i,j);
        }
    }

    std::cout << IO::done;
    return B;
}

auto solve_hydrogen_like(Atom::AtomicSystem &atom, Wavefunction &psi, SplineBasis &basis) -> void
{
    auto Hamiltonian = construct_hamiltonian_matrix(atom, psi, basis);
    auto BMatrix = construct_spline_matrix(basis);

    //SquareMatrix<double> eigenvectors(Hamiltonian.get_size());
    //std::vector<double> energies(Hamiltonian.get_size());

    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, BMatrix);

    auto N = eigenvectors.get_size();

    int k = 10;

    std::vector<double> predicted(k);
    for (int i = 1; i < k; ++i) 
        predicted.at(i-1) = atom.Z * atom.Z / (2.0 * i * i);

    
    std::cout << "---- Energies ----" << '\n';
    std::cout << " n | En    | predicted\n"; 
    for (int i=1; i < k; ++i)
    {
        printf("%2d | %1.3f | %1.3f\n", i, energies.at(i-1), predicted.at(i-1));
    }

    psi.expansion_coeffs.resize(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        psi.expansion_coeffs.at(i) = (eigenvectors(psi.n-1, i)); // check this
    }

    psi.P.resize(basis.r_grid.size());
    psi.amplitude.resize(basis.r_grid.size());

    #pragma omp parallel for
    for (std::size_t i = 0; i < basis.r_grid.size(); ++i)
    {
        for(int j = 0; j < basis.num_spl; ++j)
        {
            psi.P.at(i) += basis.bspl.at(j).at(i) * psi.expansion_coeffs.at(j);
        }
        psi.amplitude.at(i) = psi.P.at(i) / basis.r_grid.at(i);
    }

    std::cout << "\n\u222Bdr|P(r)|\u00B2 = " << trapz_linear(basis.dr, psi.P * psi.P) << ".\n";

}