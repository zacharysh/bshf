#include "solve.hpp"


auto construct_hamiltonian(const Atom &atom, const int l_state) -> SquareMatrix<double>
{
    IO::log_params(LogType::info, "Constructing Hamiltonian matrix", {{"l", l_state}});

    SquareMatrix<double> H = atom.kinetic;

    // Electrostatic nuclear potential plus the centrifugal term.  
    std::vector<double> potential = atom.nuclear_potential.values + 0.5 * l_state * (l_state + 1.0) * atom.get_range_inv() * atom.get_range_inv();

    if ( !atom.interaction_potential.values.empty() )
    {
        potential += atom.interaction_potential.values;
    }

    // We only fill the bottom half since DSYGV_ anticipates a guaranteed symmetric-definite matrix.
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) += simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j));
        }
    }

    IO::done();
    return H;
}

auto solve_schrodinger_state(Atom &atom, const int n, const int l_state) -> Electron
{
    IO::log_params(LogType::info, "H-like radial equation", {{"n", n}, {"l", l_state}}, 1);

    auto Hamiltonian = construct_hamiltonian(atom, l_state);

    //auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    Electron psi = {n, l_state, energies.at(n - l_state - 1), eigenvectors.get_row(n - l_state - 1), atom.basis};

    IO::done(-1);
    return psi;
}

auto solve_schrodinger(Atom &atom, const int l_state) -> void
{
    IO::log_params(LogType::info, "H-like radial equation", {{"l", l_state}}, 1);

    auto Hamiltonian = construct_hamiltonian(atom, l_state);

    //auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);
    
    for(int n = l_state + 1; n <= atom.n_max; ++n)
    {
        Electron psi {n, l_state, energies.at(n - l_state - 1), eigenvectors.get_row(n - l_state - 1), atom.basis};
        atom.electrons.reserve(atom.electrons.size() + sizeof(Electron));
        atom.electrons.push_back(psi);
    }

    IO::done(-1);
}

auto solve_atom(Atom &atom, int l_max) -> void
{
    IO::log("Constructing atomic spectrum", 1);

    for (int l = 0; l <= atom.l_max; ++l)
    {
        solve_schrodinger(atom, l);
    }

    // Better to sort by ascending order of energies.
    std::sort(atom.electrons.begin(), atom.electrons.end(),
                        [] (const Electron &a, const Electron &b) { return a.energy < b.energy; });

    IO::done(-1);
}

auto solve_excited_valence(Atom &atom, const int n_state, const int l_state) -> void
{
    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    IO::log_params(LogType::info, "Constructing excited valence state", {{"n", n_state}, {"l", l_state}}, 1);
    
    atom.valence().filled = false;

    if(atom.interaction_potential.type == Potential::Type::HartreeFock)
        atom.electrons.push_back(std::move(HartreeFock::solve_full_schrodinger_state(atom, n_state, l_state)));
    else
        atom.electrons.push_back(std::move(solve_schrodinger_state(atom, n_state, l_state)));

    IO::done(-1);
}

auto greens_perturbation(const Atom &atom, Electron &psi) -> void
{
    // We want to know that perturbation theory has been applied so that we can print the correction.
    psi.has_correction = true;

    // Compute the expectation value of the Greens interaction, 
    // i.e. <V_gr> = \int_0^\infty |P|^2 V_gr dr.
    IO::log("Computing Greens expectation value");
    auto greens_correction = simpson_linear(atom.basis.r_grid.dr, psi.P * psi.P * atom.interaction_potential.values);
    IO::done();

    // Now compute the e-e perturbative contribution. See Eq. 19 of task sheet.
    // We neglect the exchange term.
    IO::log("Computing direct interaction");
    auto y0ab = 2.0 * ykab(0, atom.core().P, atom.core().P, atom.basis.r_grid);
    IO::done();

    auto direct_correction = simpson_linear(atom.basis.r_grid.dr, psi.P * psi.P * y0ab);

    // Store the energy correction.
    psi.energy_correction = direct_correction - greens_correction;

    //IO::msg::done();
}