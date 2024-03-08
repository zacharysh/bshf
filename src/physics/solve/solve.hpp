#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>

#include <algorithm> //std::sort

#include <cassert> // assert

#include "../../math/calculateYK.hpp"

#include "../../math/matrix.hpp"
#include "../../math/basis/spline_basis.hpp"
#include "../../math/basis/grid.hpp"
#include "../../math/math_other.hpp"

#include "../../IO/io.hpp"

#include "../atom.hpp"
#include "../wavefunction.hpp"
#include "../potential.hpp"


auto construct_hamiltonian_matrix(const Atom &system, const int l_number) -> SquareMatrix<double>;

auto solve_schrodinger_state(Atom &atom, const int n_number, const int l_number) -> Electron;

auto solve_schrodinger(Atom &atom, const int l_number, const bool consider_atom_Z) -> void;

auto solve_atom(Atom &atom) -> void;

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void;

auto greens_perturbation(const Atom &atom, Electron &psi) -> void;


#endif