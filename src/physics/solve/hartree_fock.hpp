#ifndef HARTREE_FOCK_HPP_
#define HARTREE_FOCK_HPP_

#include <vector>

#include "../../math/matrix.hpp"

#include "solve.hpp"

#include "../potential.hpp"
#include "../atom.hpp"
#include "../wavefunction.hpp"

namespace HartreeFock
{

auto construct_full_hamiltonian(const Atom &atom, int l_number) -> SquareMatrix<double>;

auto solve_full_schrodinger_state(const Atom &atom, int n, int l_number) -> Electron;
auto solve_full_schrodinger(Atom &atom, int l_number)  -> void;

auto procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>;

auto solve_self_consistent(Atom &atom) -> void;

auto solve(Atom &atom) -> void;

auto solve_excited_valence(Atom &atom, int n_number, int l_number) -> void;
auto solve_full_excited_valence(Atom &atom, int n_number, int l_number) -> void;

};

#endif