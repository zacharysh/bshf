#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>
#include <algorithm> //std::sort

#include "../../io.hpp"

#include "../../math/matrix.hpp"
#include "../../math/basis/spline_basis.hpp"
#include "../../math/math_other.hpp"

#include "../atom.hpp"
#include "../wavefunction.hpp"
#include "../potential.hpp"


auto construct_hamiltonian(const Atom &atom, const int l_state) -> SquareMatrix<double>;

auto solve_schrodinger_state(Atom &atom, const int n_state, const int l_state) -> Electron;

auto solve_schrodinger(Atom &atom, const int l_number) -> void;

auto solve_atom(Atom &atom, const int l_max = 0) -> void;

auto solve_excited_valence(Atom &atom, const int n, const int l) -> void;

auto greens_perturbation(const Atom &atom, Electron &psi) -> void;


#endif