#ifndef HARTREE_FOCK_HPP_
#define HARTREE_FOCK_HPP_

#include <vector>

#include "../../io.hpp"

#include "../../math/matrix.hpp"

#include "solve.hpp"

#include "../atom.hpp"
#include "../potential.hpp"
#include "../wavefunction.hpp"

namespace HartreeFock
{

auto construct_hamiltonian(const Atom &atom, const int l_state, const int l_max) -> SquareMatrix<double>;

auto solve_full_schrodinger_state(const Atom &atom, const int n, const int l_number) -> Electron;

auto solve_full_schrodinger(Atom &atom)  -> void;

auto hartree_procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>;

auto solve_self_consistent(Atom &atom) -> void;

auto solve(Atom &atom) -> void;

//auto solve_excited_valence(Atom &atom, int n_number, int l_number) -> void;
//auto solve_full_excited_valence(Atom &atom, int n_number, int l_number) -> void;

};

#endif