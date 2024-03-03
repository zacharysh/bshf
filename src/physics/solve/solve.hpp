#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>

#include <algorithm> //std::sort

#include "../../math/matrix.hpp"
#include "../../math/basis/spline_basis.hpp"
#include "../../math/basis/grid.hpp"
#include "../../math/math_other.hpp"

#include "../../IO/io.hpp"

#include "../atom.hpp"
#include "../wavefunction.hpp"


auto construct_hamiltonian_matrix(const Atom &system, int l_number) -> SquareMatrix<double>;

auto solve_schrodinger(Atom &atom, int l_number, int n_number) -> void;
auto solve_schrodinger(Atom &atom, int l_number, bool consider_atom_Z) -> void;

auto solve_atom(Atom &atom) -> void;

auto solve_hartree_fock(Atom &atom) -> void;


#endif