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

    // We only fill the bottom half since DSYGV_ anticipates a symmetric-definite matrix.
    #pragma omp parallel for
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) = trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl_derivative.at(i) * atom.basis.bspl_derivative.at(j)) / 2.0
                    + trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    + trapz_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * centrifugal * atom.basis.bspl.at(j));
        }
    }

    IO::msg::done();
    return H;
}

// Duplicate.......
auto solve_schrodinger(Atom &atom, int l_number, int n) -> void
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"n", n}, {"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, l_number);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    
    // Check n >= l + 1.
        
    Electron psi{n, l_number, 0, energies.at(n - 1), eigenvectors.get_row(n - 1), atom.basis};
    atom.electrons.reserve(atom.electrons.size() + sizeof(Electron));
    atom.electrons.push_back(psi);
    
    IO::msg::done(true);
}

auto solve_schrodinger(Atom &atom, int l_number, bool consider_atom_Z) -> void
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, l_number);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    
    for(int n = l_number + 1; n <= (l_number + 1) + 2; ++n)
    {
        if((2 * static_cast<int>(atom.electrons.size()) >= atom.Z) && consider_atom_Z)
            break;
        
        Electron psi{n, l_number, 0, energies.at(n - 1), eigenvectors.get_row(n - 1), atom.basis};
        atom.electrons.reserve(atom.electrons.size() + sizeof(Electron));
        atom.electrons.push_back(psi);
    }
    IO::msg::done(true);

}

auto solve_atom(Atom &atom) -> void
{
    IO::msg::action("Constructing", "atomic spectra", true);
    // How to do this?
    for (int l = 0; l < 3; ++l)
    {
        solve_schrodinger(atom, l, true);
    }

    std::sort(atom.electrons.begin(), atom.electrons.end(),
                        [](const Electron &a, const Electron &b) { return a.energy < b.energy; });
    
    atom.electrons.resize((atom.Z + 1) / 2);

    IO::msg::done();
}

auto solve_hartree_fock(Atom &atom) -> void
{
    IO::msg::action("Solve", "Hartree-Fock Method", true);

    IO::msg::done();
}