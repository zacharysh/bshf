#include "hartree_fock.hpp"


namespace HartreeFock
{

// Requires states to be done already.
auto solve_full_schrodinger(Atom &atom, int l_number) -> void
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"l", l_number}}, true);
    
    auto Hamiltonian = construct_full_hamiltonian(atom, l_number, atom.electrons.back().l);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    
    for(std::size_t i = 0; i < atom.electrons.size(); ++i)
    {        
        if (atom.electrons.at(i).l == l_number)
        {
            auto n = atom.electrons.at(i).n;

            atom.electrons.at(i) = {n, l_number, 0, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis};

        }
    }
    IO::msg::done(true);

}

auto solve_full_schrodinger_state(const Atom &atom, int n, int l_number) -> Electron
{
    IO::msg::action<int>("Solving", "H-like radial equation", {{"n", n}, {"l", l_number}}, true);

    auto Hamiltonian = construct_full_hamiltonian(atom, l_number, atom.electrons.back().l);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    auto psi = Electron(n, l_number, 0, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis);

    IO::msg::done(true);
    return psi;
}

auto construct_full_hamiltonian(const Atom &atom, const int l_state, const int l_max) -> SquareMatrix<double>
{
    IO::msg::construct<int>("Hamiltonian matrix with exchange term", {{"l", l_state}});

    // generate centrifugal term
    std::vector<double> centrifugal(atom.basis.grid_size());
    std::transform(atom.get_range_inv().begin(), atom.get_range_inv().end(), centrifugal.begin(),
            [l_state] (auto r_inv) { return 0.5 * l_state * (l_state + 1.0) * r_inv * r_inv; });

    SquareMatrix<double> H(atom.basis.num_spl);

    auto potential = atom.nuclear_potential.values + centrifugal + atom.interaction_potential.values;

    // We only fill the bottom half since DSYGV_ anticipates a guaranteed symmetric-definite matrix.

    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            auto lambda = (l_max == 0) ? 0.5 : 1.0/6.0;
            std::vector<double> exchange_j = -2.0 * lambda * YK::ykab(l_max, atom.electrons.front().P, atom.basis.bspl.at(j), atom.get_range()) * atom.electrons.front().P;
            
            H(i, j) = 0.5 * simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl_derivative.at(i) * atom.basis.bspl_derivative.at(j))
                    + simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    + simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * exchange_j);
        }
    }

    IO::msg::done();
    return SquareMatrix<double>(H);
}

auto procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>
{
    std::vector<double> core_energies;
    core_energies.reserve(20); // Typically takes this long.

    auto iteration = 0;
    auto rel_energy_diff = 1.0; // Initialise?

    core_energies.push_back(atom.electrons.front().energy);

    const auto energy_threshold = 1e-6;

    

    while (rel_energy_diff >= energy_threshold)
    {
        IO::msg::action<double>("Performing", "Hartree iteration", {{"iteration", iteration}});

        atom.interaction_potential.values = 2.0 * YK::ykab(0, atom.electrons.front().P, atom.electrons.front().P, atom.get_range());

        // Be quiet and save screen space.
        IO::msg::verbose = false;

        if(full_hamiltonian)
        {
                //atom.electrons.front() = solve_full_schrodinger(atom, iter.l);
            for (auto iter : atom.electrons){
                solve_full_schrodinger(atom, iter.l);}
        }
        else
        {
            atom.electrons.front() = solve_schrodinger_state(atom, 1, 0);
        }

        // Unnecessary copying, but I think makes the algorithm easier to read.
        auto prev_core_energy = core_energies.back();
        core_energies.push_back(atom.electrons.front().energy);

        rel_energy_diff = abs(abs(core_energies.back()) - abs(prev_core_energy)) / abs(prev_core_energy);

        ++iteration;

        
        // We can be noisy again.
        IO::msg::verbose = true;

        IO::msg::done();

        // Don't want to overflow.
        if(iteration >= 20)
            core_energies.reserve(core_energies.size() + sizeof(double));
    }

    return core_energies;
}

auto solve_self_consistent(Atom &atom) -> void
{
    IO::msg::action("Solving", "atom with Self-Consistent-Hartree-Fock interaction", true);
    
    // First, solve Schrodinger equation with Greens.
    atom.interaction_potential = { Potential(atom.Z, Potential::Type::Greens, atom.basis.r_grid) };
    solve_atom(atom);

    std::vector<double> initial_energies(atom.get_energies());

    // Book-keeping
    atom.interaction_potential.type = Potential::Type::HF_Direct;

    // Compute Hartree convergence algorithm.
    auto core_energy_timeseries = procedure(atom, false);

    atom.electrons = {};

    // Finally, solve the full atom with the new direct potential.
    solve_atom(atom);

    IO::msg::done(true);
}

auto solve(Atom &atom) -> void
{
    // Only valid for Lithium at the moment (i.e. only 1s and 2s states are stored).
    assert(atom.Z == 3 && "Only works for lithium at the moment.");

    IO::msg::action("Solving", "atom with Hartree-Fock method", true);
    
    // Start with self-consistent solution.
    solve_self_consistent(atom);
    std::vector<double> initial_energies(atom.get_energies());

    // Compute Hartree-Fock convergence algorithm.
    auto core_energy_timeseries = procedure(atom, true);


    IO::msg::done(true);
}

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void
{
    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    IO::msg::construct<int>("excited valence state", {{"n", n}, {"l", l}}, true);
    
    // Assuming the valence is the back, which it should be.
    atom.electrons.back().filled = false;
    atom.electrons.push_back(solve_full_schrodinger_state(atom, n, l));

    IO::msg::done(true);
}

auto solve_full_excited_valence(Atom &atom, const int n, const int l) -> void
{
    // Only valid for Lithium at the moment (i.e. only 1s and 2s states are stored).
    assert((n == 2 && l == 1) && "Only works for lithium at the moment.");

    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    auto empty_shell = std::move(atom.electrons.back());

    // Assuming the valence is the back, which it should be.
    atom.electrons.back() = solve_full_schrodinger_state(atom, n, l);
    
    // Rerun Hartree-Fock with 2p instead of 1s.
    auto core_energy_timeseries = procedure(atom, true);

    auto valence_shell = std::move(atom.electrons.back());
    

    atom.electrons.back() = std::move(empty_shell);
    atom.electrons.back().filled = false;

    atom.electrons.push_back(std::move(valence_shell));

    IO::msg::done(true);
}

}