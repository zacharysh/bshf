#include "hartree_fock.hpp"


namespace HartreeFock
{

auto construct_hamiltonian(const Atom &atom, const int l_state, const int l_max) -> SquareMatrix<double>
{
    IO::log_params(LogType::info, "Hamiltonian matrix with exchange term", {{"l", l_state}});

    SquareMatrix<double> H {atom.kinetic};

    // Electrostatic nuclear potential plus the centrifugal term.
    std::vector<double> potential = {
        atom.nuclear_potential.values + atom.interaction_potential.values + 0.5 * l_state * (l_state + 1.0) * atom.get_range_inv() * atom.get_range_inv()
        };

    auto lambda = (l_max == 0) ? 1.0/2.0 : 1.0/6.0;

    // We only fill the bottom half since dsygv_ anticipates a symmetric-definite matrix.
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        // Since H is Hermitian, we can save significant time by computing the action of V_ex b_i N times instead of V_ex b_j N^2 times.
        std::vector<double> exchange = -2.0 * lambda * ykab(l_max, atom.core().P, atom.basis.bspl.at(i), atom.basis.r_grid) * atom.core().P;
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) += simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    +  simpson_linear(atom.basis.r_grid.dr,                         exchange  * atom.basis.bspl.at(j));
        }
    }

    IO::done();
    return SquareMatrix<double>(H);
}

// Requires states to be done already.
auto solve_full_schrodinger(Atom &atom, int l_number) -> void
{
    IO::log_params(LogType::info, "Solving H-like radial equation", {{"l", l_number}}, 1);
    
    auto Hamiltonian = construct_hamiltonian(atom, l_number, atom.electrons.back().l);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    
    for(std::size_t i = 0; i < atom.electrons.size(); ++i)
    {        
        if (atom.electrons.at(i).l == l_number)
        {
            auto n = atom.electrons.at(i).n;

            atom.electrons.at(i) = {n, l_number, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis};

        }
    }
    
    IO::done(-1);
}

// Requires states to be done already.
auto solve_full_schrodinger(Atom &atom) -> void
{    
    for(Electron &psi : atom.electrons)
    {        
        IO::log_params(LogType::info, "Solving H-like radial equation", {{"l", psi.l}}, 1);
        auto Hamiltonian = construct_hamiltonian(atom, psi.l, atom.valence().l);
        auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);

        psi = {psi.n, psi.l, energies.at(psi.n - psi.l - 1), eigenvectors.get_row(psi.n - psi.l - 1), atom.basis};
        IO::done(-1);
    }
}


auto solve_full_schrodinger_state(const Atom &atom, int n, int l_number) -> Electron
{
    IO::log_params(LogType::info, "Solving H-like radial equation", {{"n", n}, {"l", l_number}}, 1);

    auto Hamiltonian = construct_hamiltonian(atom, l_number, atom.valence().l);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(Hamiltonian), atom.basis.Bmatrix);
    Electron psi {n, l_number, energies.at(n - l_number - 1), eigenvectors.get_row(n - l_number - 1), atom.basis};

    IO::done(-1);
    return psi;
}





auto hartree_procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>
{
    // Typically takes this many iterations.
    std::vector<double> core_energies;
    core_energies.reserve(20);

    auto iteration = 0;
    auto rel_energy_diff = 1.0; // Initialise?

    core_energies.push_back(atom.electrons.front().energy);

    // Iterate until relative core energy change is negligible.
    while (rel_energy_diff >= 1e-6)
    {
        IO::log_params(LogType::info, "Performing Hartree procedure", {{"iteration", iteration}});

        atom.interaction_potential.values = {2.0 * ykab(0, atom.core().P, atom.core().P, atom.basis.r_grid)};

        // Be quiet and save screen space.
        IO::verbose = false;

        if(full_hamiltonian)
        {
            solve_full_schrodinger(atom);
        }
        else
        {
            atom.core() = solve_schrodinger_state(atom, 1, 0);
        }

        // Unnecessary copying, but I think makes the algorithm easier to read.
        auto prev_core_energy = core_energies.back();
        core_energies.push_back(atom.core().energy);

        rel_energy_diff = abs(abs(core_energies.back()) - abs(prev_core_energy)) / abs(prev_core_energy);

        ++iteration;

        
        // We can be noisy again.
        IO::verbose = true;
        

        IO::done();

        // Don't want to overflow.
        if(iteration >= 20)
            core_energies.reserve(core_energies.size() + sizeof(double));
    }
    

    return core_energies;
}

auto solve_self_consistent(Atom &atom) -> void
{
    IO::log("Solving atom with mean-field Hartree approximation", 1);
    
    // First, solve Schrodinger equation with Greens.
    atom.interaction_potential = {atom.Z, Potential::Type::Greens, atom.basis.r_grid};
    solve_atom(atom);

    std::vector<double> initial_energies(atom.get_energies());

    // Book-keeping
    atom.interaction_potential.type = Potential::Type::HF_Direct;

    // Compute Hartree convergence algorithm.
    auto core_energy_timeseries = hartree_procedure(atom, false);

    atom.electrons = {};

    // Finally, solve the full atom with the new direct potential.
    solve_atom(atom);

    IO::done(-1);
}

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void
{
    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    IO::log_params(LogType::info, "excited valence state", {{"n", n}, {"l", l}}, 1);
    
    // Assuming the valence is the back, which it should be.
    atom.electrons.back().filled = false;
    atom.electrons.push_back(solve_full_schrodinger_state(atom, n, l));

    IO::done(-1);
}

auto solve(Atom &atom) -> void
{
    // Only valid for Lithium at the moment (i.e. only 1s and 2s states are stored).
    assert(atom.Z == 3 && "Only works for lithium at the moment.");

    IO::log("Solving atom with Hartree-Fock method", 1);
    
    // Start with self-consistent solution.
    solve_self_consistent(atom);
    std::vector<double> initial_energies(atom.get_energies());

    // Compute Hartree-Fock convergence algorithm.
    auto core_energy_timeseries = hartree_procedure(atom, true);


    IO::done(-1);
}

auto solve_full_excited_valence(Atom &atom, const int n, const int l) -> void
{
    IO::log_params(LogType::info, "Computing excited valence state", {{"n", n}, {"l", l}}, 1);

    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

    auto empty_shell = std::move(atom.valence());

    // Assuming the valence is the back, which it should be.
    atom.electrons.back() = solve_full_schrodinger_state(atom, n, l);
    
    // Rerun Hartree-Fock with 2p instead of 1s.
    auto core_energy_timeseries = hartree_procedure(atom, true);

    auto valence_shell = std::move(atom.valence());
    
    atom.valence() = std::move(empty_shell);
    atom.valence().filled = false;

    atom.electrons.push_back(std::move(valence_shell));

    IO::done(-1);
}

}