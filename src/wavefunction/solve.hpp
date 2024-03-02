#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>

#include "../math/matrix.hpp"
#include "../math/SplineBasis.hpp"
#include "../math/math_other.hpp"

#include "../IO/io.hpp"

#include "atom.hpp"
#include "wavefunction.hpp"


auto construct_hamiltonian_matrix(Atom::AtomicSystem system, Wavefunction psi, const SplineBasis &basis) -> SquareMatrix<double>;

auto construct_spline_matrix(const SplineBasis &basis) -> SquareMatrix<double>;

auto solve_hydrogen_like(Atom::AtomicSystem &atom, Wavefunction &psi, SplineBasis &basis) -> void;


#endif