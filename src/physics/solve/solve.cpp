#include "solve.hpp"


auto construct_hamiltonian_matrix(const Atom &atom, int l_number) -> SquareMatrix<double>
{
    IO::msg::construct<int>("Hamiltonian matrix", {{"l", l_number}});

    // generate centrifugal term
    std::vector<double> centrifugal(atom.basis.grid_size());

    for (std::size_t i = 0; i < atom.basis.grid_size(); ++i) 
        centrifugal.at(i) =  l_number * (l_number + 1)  / (2.0 * atom.basis.r_grid.at(i) * atom.basis.r_grid.at(i));

    SquareMatrix<double> H(atom.basis.num_spl);

    std::vector<double> potential = atom.nuclear_potential.get_values();

    // AAAHHHHH! Fix me!
    if ( !atom.interaction_potential.get_values().empty() )
    {
        potential = potential + atom.interaction_potential.get_values(); // Fix operators.
    }
    else
        std::cout << " (no interaction)";

    // We only fill the bottom half since DSYGV_ anticipates a guaranteed symmetric-definite matrix.
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

auto solve_schrodinger_state(const Atom &atom, int n, int l_number) -> Electron
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"n", n}, {"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, l_number);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    auto psi = Electron(n, l_number, 0, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis);

    //IO::msg::print_values<double>({{"energy", psi.energy}});

    IO::msg::done(true);
    return psi;
}


// Delete?
auto solve_schrodinger(Atom &atom, int l_number, bool consider_atom_Z) -> void
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"l", l_number}}, true);

    auto Hamiltonian = construct_hamiltonian_matrix(atom, l_number);

    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    
    
    for(int n = l_number + 1; n <= (l_number + 1) + 2; ++n)
    {
        if((2 * static_cast<int>(atom.electrons.size()) >= atom.Z) && consider_atom_Z)
            break;
        
        Electron psi{n, l_number, 0, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis};
        atom.electrons.reserve(atom.electrons.size() + sizeof(Electron));
        atom.electrons.push_back(psi);
    }
    IO::msg::done(true);

}

auto solve_atom(Atom &atom) -> void
{
    IO::msg::action("Constructing", "atomic spectra", true);

    // FIX ME!
    // How to do this?
    for (int l = 0; l <= 1; ++l)
    {
        solve_schrodinger(atom, l, true);
    }

    std::sort(atom.electrons.begin(), atom.electrons.end(),
                        [](const Electron &a, const Electron &b) { return a.energy < b.energy; });
    
    // Only valid for Lithium at the moment (i.e. only 1s and 2s states are stored.).
    atom.electrons.resize(2);

    IO::msg::done(true);
}

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void
{
    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    IO::msg::construct<int>("excited valence state", {{"n", n}, {"l", l}}, true);
    
    // Assuming the valence is the back, which it should be.
    atom.electrons.back().filled = false;
    atom.electrons.push_back(solve_schrodinger_state(atom, n, l));

    IO::msg::done(true);
}

auto greens_perturbation(const Atom &atom, Electron &psi) -> void
{
    // We want to know that perturbation theory has been applied so that we can print the correction.
    // Perhaps this is redundant since this is stored by the atomic interaction_potential parameter?
    psi.has_correction = true;

    // Compute the expectation value of the Greens interaction, 
    // i.e. <V_gr> = \int_0^\infty |P|^2 V_gr dr.
    IO::msg::action("Solving", "Greens expectation");
    auto greens_correction = trapz_linear(atom.basis.r_grid.dr, psi.P * psi.P * atom.interaction_potential.get_values());
    IO::msg::done();

    // Now compute the e-e perturbative contribution. See Eq. 19 of task sheet.
    // We neglect the exchange term.
    // Need to generalise to all?
    assert(atom.Z == 3 && "Only working for Li, currently.");

    IO::msg::action("Solving", "direct contribution");
    auto y0ab = 2.0 * YK::ykab(0, atom.electrons.front().P, atom.electrons.front().P, atom.basis.r_grid.range);
    IO::msg::done();

    auto direct_correction = trapz_linear(atom.basis.r_grid.dr, psi.P * psi.P * y0ab);

    // Store the energy correction.
    psi.energy_correction = direct_correction - greens_correction;

    //IO::msg::done();
}