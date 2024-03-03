#include "solve.hpp"


auto construct_hamiltonian_matrix(const Atom &atom, int l_number) -> SquareMatrix<double>
{
    IO::msg::construct<int>("Hamiltonian matrix", {{"l", l_number}});

    // generate centrifugal term
    std::vector<double> centrifugal(atom.basis.grid_size());

    for (std::size_t i = 0; i < atom.basis.grid_size(); ++i) 
        centrifugal.at(i) =  l_number * (l_number + 1)  / (2.0 * atom.basis.r_grid.at(i) * atom.basis.r_grid.at(i));

    auto potential(atom.nuclear_potential.get_values());

    // Check if interaction is present.
    if(!atom.interaction_potential.get_values().empty())
        potential = potential + atom.interaction_potential.get_values(); // Fix operators


    SquareMatrix<double> H(atom.basis.num_spl);

    #pragma omp parallel for
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) = trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl_derivative.at(i) * atom.basis.bspl_derivative.at(j)) / 2.0
                    + trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    + trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * centrifugal * atom.basis.bspl.at(j));

            // Only need to calculate bottom half.
            //H(j,i) = H(i,j);
        }
    }

    IO::msg::done();
    return H;
}

auto solve_schrodinger(Atom &atom, int l_number) -> void
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, l_number);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    
    for(int n = l_number + 1; n <= (l_number + 1) + 2; ++n)
    {
        //if(2 * static_cast<int>(atom.electrons.size()) >= atom.Z)
        //    break;
        
        Electron psi{n, l_number, 0, energies.at(n - 1), eigenvectors.get_row(n - 1), atom.basis};
        atom.electrons.reserve(atom.electrons.size() + sizeof(Electron));
        atom.electrons.push_back(psi);
    }
    IO::msg::done(true);

}

auto solve_atom(Atom &atom) -> void
{
    int l = 0;
    while(2 * static_cast<int>(atom.electrons.size()) <= atom.Z + 1)
    {
        solve_schrodinger(atom, l);
        ++l;
    }
}