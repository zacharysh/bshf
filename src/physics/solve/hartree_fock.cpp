#include "hartree_fock.hpp"


namespace HartreeFock
{

auto construct_hamiltonian(const Atom &atom, const int l_state, const int n_max) -> SquareMatrix<double>
{
    IO::log_params(LogType::info, "Hamiltonian matrix with exchange term", {{"l", l_state}});

    SquareMatrix<double> H = atom.kinetic;

    // Electrostatic nuclear potential AND direct HF term plus the centrifugal term.
    std::vector<double> potential = 
        atom.nuclear_potential.values + atom.interaction_potential.values + 0.5 * l_state * (l_state + 1.0) * atom.get_range_inv() * atom.get_range_inv();

    auto lambda2 = (l_state == 0) ? -1.0 : -1.0/3.0;

    // We only fill the bottom half since dsygv_ anticipates a symmetric-definite matrix.
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        // Since H is Hermitian, we can save significant time by computing the action of V_ex b_i N times instead of V_ex b_j N^2 times.
        std::vector<double> exchange = ykab(l_state, atom.core().P, atom.basis.bspl.at(i), atom.basis.r_grid) * atom.core().P;
        
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) += simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    +  simpson_linear(atom.basis.r_grid.dr,                         exchange  * atom.basis.bspl.at(j)) * lambda2;
        }
    }

    IO::done();
    return SquareMatrix<double>(H);
}


// Requires states to be done already.
auto solve_full_schrodinger(Atom &atom) -> void
{    
    for(Electron &psi : atom.electrons)
    {
        IO::log_params(LogType::info, "Solving radial equation", {{"l", psi.l}}, 1);
        auto Hamiltonian = construct_hamiltonian(atom, psi.l, psi.n);//atom.valence().l);
        auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);

        psi = Electron(psi.n, psi.l, energies.at(psi.n - psi.l - 1), eigenvectors.get_row(psi.n - psi.l - 1), atom.basis);
        IO::done(-1);
    }
}

auto solve_full_schrodinger_state(const Atom &atom, int n, int l) -> Electron
{
    IO::log_params(LogType::info, "Solving radial equation", {{"n", n}, {"l", l}}, 1);

    auto Hamiltonian = construct_hamiltonian(atom, l, n);//atom.valence().l);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(Hamiltonian, atom.basis.Bmatrix);

    IO::done(-1);
    return Electron(n, l, energies.at(n - l - 1), eigenvectors.get_row(n - l - 1), atom.basis);
}

auto hartree_procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>
{
    // Typically takes this many iterations.
    std::vector<double> core_energies;
    core_energies.reserve(20);

    auto iteration = 0;
    auto rel_energy_diff = 1.0; // Initialise?

    core_energies.push_back(atom.core().energy);

    // To return to the initial verbosity, in case this is set to false initially.
    auto initial_verbosity = IO::verbose;

    // Iterate until relative core energy change is negligible.
    while (rel_energy_diff >= 1e-6)
    {
        IO::log_params(LogType::info, "Performing Hartree procedure", {{"iteration", iteration}});

        atom.interaction_potential.values = 2.0 * ykab(0, atom.core().P, atom.core().P, atom.basis.r_grid);

        // Be quiet and save screen space.
        IO::verbose = false;

        if(full_hamiltonian)
        {
            atom.core() = solve_full_schrodinger_state(atom, 1, 0);
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

        
        // We can be noisy again (if we want).
        IO::verbose = initial_verbosity;

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

auto solve(Atom &atom) -> void
{
    IO::log("Solving atom with Hartree-Fock method", 1);
    
    // Start with self-consistent solution.
    solve_self_consistent(atom);
    std::vector<double> initial_energies(atom.get_energies());

    // Compute Hartree-Fock convergence algorithm.
    auto core_energy_timeseries = hartree_procedure(atom, true);

    solve_full_schrodinger(atom);

    atom.interaction_potential.type = Potential::Type::HartreeFock;

    IO::done(-1);
}

}