#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>    //std::vector
#include <algorithm> //std::sort

#include "../../io.hpp"

#include "../../math/matrix.hpp"
#include "../../math/basis/spline_basis.hpp"
#include "../../math/math_other.hpp"

#include "../atom.hpp"
#include "../electron.hpp"
#include "../potential.hpp"
#include "./hartree_fock.hpp"


auto construct_hamiltonian(const Atom &atom, const int l_state) -> SquareMatrix<double>;
auto schrodinger(Atom &atom, const int l_number) -> void;
auto solve_schrodinger_state(const Atom &atom, const int n_state, const int l_state) -> Electron;

auto greens_perturbation(const Atom &atom, Electron &psi) -> void;
auto generate_atom(Atom &atom, bool compute_perturbation = true) -> void;
auto solve_excited_valence(Atom &atom, const int n, const int l) -> void;

#endif