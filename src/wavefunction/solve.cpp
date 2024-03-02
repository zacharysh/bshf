#include "solve.hpp"


auto construct_hamiltonian_matrix(const Atom::AtomicSystem &system, const SplineBasis &basis, int l_number) -> SquareMatrix<double>
{
    IO::msg::construct<int>("Hamiltonian matrix", {{"l", l_number}});

    
    // fix me up
    // generate centrifugal term
    std::vector<double> centrifugal(basis.r_grid.size());

    for (std::size_t i = 0; i < basis.r_grid.size(); ++i) 
        centrifugal.at(i) =  l_number * (l_number + 1)  / (2.0 * basis.r_grid.at(i) * basis.r_grid.at(i));


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
    IO::msg::done();
    return H;
}



auto solve_hydrogen_like(Atom::AtomicSystem &atom, SplineBasis &basis, int l_number, int max_n_number) -> std::vector<Electron>
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, basis, l_number);

    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, basis.Bmatrix);

    auto N_b = eigenvectors.get_size();

    std::vector<Electron> solutions {};

    solutions.reserve(N_b);

    for (int n = 1; n < max_n_number; ++n)
    {
        Electron psi{n, l_number, 0, energies.at(n-1), eigenvectors.get_row(n - 1), basis};
        solutions.push_back(psi);
    }

    IO::msg::done(true);
    return solutions;
}