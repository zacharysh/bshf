#include "solve.hpp"

auto construct_hamiltonian(const Atom &atom, const int l) -> SquareMatrix<double>
{
    IO::log_params(LogType::info, "Constructing Hamiltonian matrix", {{"l", l}});
    SquareMatrix<double> H {atom.kinetic};

    // Electrostatic nuclear potential and interaction potential (if present) plus the centrifugal term.  
    std::vector<double> potential {atom.nuclear_potential.values + atom.interaction_potential.values + 0.5 * l * (l + 1.0) * atom.get_range_inv() * atom.get_range_inv()};

    // We only fill the bottom half since dsygv_ anticipates a guaranteed symmetric-definite matrix.
    #pragma omp parallel for
    for(int i = 0; i < atom.basis.num_spl; ++i)
        for(int j = 0; j <= i; ++j)
            H(i, j) += simpson_linear(atom.basis.r_grid.dr, atom.basis.bspl.at(i) * potential * atom.basis.bspl.at(j));

    IO::done();
    return SquareMatrix<double>(H);
}

auto solve_schrodinger_state(const Atom &atom, const int n, const int l) -> Electron
{
    IO::log_params(LogType::info, "H-like radial equation", {{"n", n}, {"l", l}}, 1);

    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(construct_hamiltonian(atom, l)), atom.basis.Bmatrix);

    IO::done(-1);
    return Electron(n, l, energies.at(n - l - 1), eigenvectors.get_row(n - l - 1), atom.basis);
}

auto schrodinger(Atom &atom, const int l) -> void
{
    IO::log_params(LogType::info, "H-like radial equation", {{"l", l}}, 1);

    auto [eigenvectors, energies] = MatrixTools::solve_eigen_system(std::move(construct_hamiltonian(atom, l)), atom.basis.Bmatrix);
    
    for(auto n = l + 1; n <= atom.n_max; ++n)
        atom.electrons.push_back(Electron(n, l, energies.at(n - l - 1), eigenvectors.get_row(n - l - 1), atom.basis));
    
    IO::done(-1);
}

auto greens_perturbation(const Atom &atom, Electron &psi) -> void
{
    IO::log("Computing Perturbative correction.", 1);
    // We want to know that perturbation theory has been applied so that we can print the correction.
    psi.has_correction = true;

    // Compute the expectation value of the Greens interaction, 
    // i.e. <V_gr> = \int_0^\infty |P|^2 V_gr dr.
    IO::log("Computing Greens expectation value");
    auto greens_correction = simpson_linear(atom.basis.r_grid.dr, psi.P * psi.P * atom.interaction_potential.values);
    IO::done();

    // Now we compute the e-e perturbative contribution. See Eq. 19 of task sheet.
    // Of course, we neglect the exchange term here.
    IO::log("Computing direct interaction");
    auto y0ab = 2.0 * ykab(0, atom.core().P, atom.core().P, atom.basis.r_grid);
    IO::done();

    auto direct_correction = simpson_linear(atom.basis.r_grid.dr, psi.P * psi.P * y0ab);

    // Store the energy correction.
    psi.energy_correction = direct_correction - greens_correction;

    IO::done(-1);
}

auto generate_atom(Atom &atom, bool compute_perturbation) -> void
{
    IO::log("Constructing atomic spectrum", 1);

    for (int l = 0; l <= atom.l_max; ++l)
        schrodinger(atom, l);
    
    // Better to sort by ascending order of energies.
    std::sort(atom.electrons.begin(), atom.electrons.end(),
                        [] (const Electron &a, const Electron &b) { return a.energy < b.energy; });


    // May not always need to compute the perturbation despite having Greens potential, such as when doing SC Hartree.
    if(atom.interaction_potential.type == Potential::Type::Greens && compute_perturbation)
        greens_perturbation(atom, atom.valence());

    IO::done(-1);
}

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void
{
    // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
    // e.g. for Lithium this is 1s2 2s1 -> 1s2 2p1).
    IO::log_params(LogType::info, "Constructing excited valence state", {{"n", n}, {"l", l}}, 1);
    
    atom.valence().filled = false;

    if(atom.interaction_potential.type == Potential::Type::HartreeFock)
        atom.electrons.push_back(std::move(HartreeFock::solve_full_schrodinger_state(atom, n, l)));
    else
    {
        atom.electrons.push_back(std::move(solve_schrodinger_state(atom, n, l)));

        if(atom.interaction_potential.type == Potential::Type::Greens)
            greens_perturbation(atom, atom.valence());
    }

    IO::done(-1);
}