#include "hartree_fock.hpp"


namespace HartreeFock
{

auto construct_hamiltonian(const Atom &atom, const int l) -> SquareMatrix<double>
{
    IO::log_params(LogType::info, "Hamiltonian matrix with exchange term", {{"l", l}});

    SquareMatrix<double> H = atom.kinetic;

    // Electrostatic nuclear potential AND direct HF term plus the centrifugal term.
    std::vector<double> potential = 
        atom.nuclear_potential.values + atom.interaction_potential.values + 0.5 * l * (l + 1.0) * atom.get_range_inv() * atom.get_range_inv();

    auto lambda2 = (l == 0) ? -1.0 : -1.0/3.0;

    // We only fill the bottom half since dsygv_ anticipates a symmetric-definite matrix.
    #pragma omp parallel for
    for(int i = 0; i < atom.basis.num_spl; ++i)
    {
        // Since H is Hermitian, we can save significant time by computing the action of V_ex b_i N times instead of V_ex b_j N^2 times.
        std::vector<double> exchange = ykab(l, atom.core().P, atom.basis.bspl.at(i), atom.basis.r_grid) * atom.core().P;
        
        for(int j = 0; j <= i; ++j)
        {
            H(i, j) += simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j))
                    +  simpson_linear(atom.basis.r_grid.dr,                         exchange  * atom.basis.bspl.at(j)) * lambda2;
        }
    }

    IO::done();
    return SquareMatrix<double>(H);
}

auto solve_full_schrodinger_state(const Atom &atom, int n, int l) -> Electron
{
    IO::log_params(LogType::info, "Solving radial equation", {{"n", n}, {"l", l}}, 1);
    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(HartreeFock::construct_hamiltonian(atom, l)), atom.basis.Bmatrix);
    IO::done(-1);
    return Electron(n, l, energies.at(n - l - 1), eigenvectors.get_row(n - l - 1), atom.basis);
}

auto hartree_procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>
{
    // Typically takes this many iterations.
    std::vector<double> core_energies;
    core_energies.reserve(20);

    auto iter_count = 1;
    auto rel_energy_diff = 1.0; // Initialise?

    core_energies.push_back(atom.core().energy);

    // Store the initial verbosity flag to restore after completion of the procedure.
    auto initial_verbosity = IO::verbose;

    // Iterate until relative core energy change is negligible.
    while (rel_energy_diff >= 1e-6)
    {
        IO::log_params(LogType::info, "Performing Hartree procedure", {{"iteration", iter_count}});

        // Set the interaction potential in the system to the Hartree direct interaction.
        atom.interaction_potential.values = 2.0 * ykab(0, atom.core().P, atom.core().P, atom.basis.r_grid);

        // Turn off Hamiltonian/solver logging output.
        IO::verbose = false;

        // Solve with the exchange term or not for Hartree-Fock or Hartree or, respectively.
        if(full_hamiltonian)
            atom.core() = solve_full_schrodinger_state(atom, 1, 0);
        else
            atom.core() = solve_schrodinger_state(atom, 1, 0);

        // Store new core energy and compute its rel. change.
        auto prev_core_energy = core_energies.back();
        core_energies.push_back(atom.core().energy);
        rel_energy_diff = abs(abs(core_energies.back()) - abs(prev_core_energy)) / abs(prev_core_energy);

        // We can be noisy again (if we want).
        IO::verbose = initial_verbosity;
        IO::done();

        ++iter_count;
    }
    
    return core_energies;
}

auto self_consistent(Atom &atom) -> void
{
    IO::log("Solving atom with mean-field Hartree model", 1);
    
    // First, solve Schrodinger equation with Greens.
    atom.interaction_potential = {atom.Z, Potential::Type::Greens, atom.basis.r_grid};
    generate_atom(atom, false);

    // Book-keeping
    atom.interaction_potential.type = Potential::Type::HF_Direct;

    // Compute Hartree convergence algorithm.
    auto core_energy_time_series = hartree_procedure(atom, false);

    // Construct vector with the iteration number to use as the horizontal axis for plotting the energy time series.
    std::vector<double> it_count (core_energy_time_series.size());
    std::iota(it_count.begin(), it_count.end(), 0);
    IO::print_to_file("self_consistent_core_energies",  {{"iteration", it_count}, {"energy", core_energy_time_series}}, 15);

    atom.electrons = {};

    // Finally, solve the full atom with the new direct potential.
    generate_atom(atom, false);

    IO::done(-1);
}

auto solve(Atom &atom) -> void
{
    IO::log("Solving system with Hartree-Fock model", 1);
    
    bool compute_full_HF = (atom.interaction_potential.type == Potential::Type::HartreeFock);

    // Start with self-consistent solution.
    self_consistent(atom);

    if( compute_full_HF == false )
        return;

    // Run Hartree-Fock convergence algorithm.
    auto core_energy_time_series = hartree_procedure(atom, true);

    std::vector<double> it_count (core_energy_time_series.size());
    std::iota(it_count.begin(), it_count.end(), 0);
    IO::print_to_file("HF_core_energies",  {{"iteration", it_count}, {"energy", core_energy_time_series}}, 15);

    for(Electron &psi : atom.electrons)
        psi = solve_full_schrodinger_state(atom, psi.n, psi.l);
    
    atom.interaction_potential.type = Potential::Type::HartreeFock;

    IO::done(-1);
}

}
