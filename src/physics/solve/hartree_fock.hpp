#ifndef HARTREE_FOCK_HPP_
#define HARTREE_FOCK_HPP_

#include <vector>
#include <numeric> // std::iota

#include "../../io.hpp"

#include "../../math/matrix.hpp"

#include "solve.hpp"

#include "../atom.hpp"
#include "../potential.hpp"
#include "../electron.hpp"

namespace HartreeFock
{
    auto construct_hamiltonian(const Atom &atom, const int l) -> SquareMatrix<double>;
    auto solve_full_schrodinger_state(const Atom &atom, const int n, const int l_number) -> Electron;
    auto hartree_procedure(Atom &atom, bool full_hamiltonian) -> std::vector<double>;
    auto self_consistent(Atom &atom) -> void;
    auto solve(Atom &atom) -> void;
};

#endif