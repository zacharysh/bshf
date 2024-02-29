#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>

#include "../math/matrix.hpp"

#include "atom.hpp"
#include "wavefunction.hpp"


#include "../math/SplineBasis.hpp"

#include "../math/math_other.hpp"




auto constructHamiltonian(Atom::AtomicSystem system, Wavefunction psi, SplineBasis &basis) -> SquareMatrix<double>;

auto constructBMatrix(SplineBasis &basis) -> SquareMatrix<double>;

auto solveHydrogenlike(Atom::AtomicSystem &atom, Wavefunction &psi, SplineBasis &basis) -> void;


#endif