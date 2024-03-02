#ifndef SOLVE_HPP_
#define SOLVE_HPP_

#include <vector>

#include "../math/matrix.hpp"
#include "../math/SplineBasis.hpp"
#include "../math/math_other.hpp"

#include "../IO/io.hpp"

#include "atom.hpp"
#include "wavefunction.hpp"


auto construct_hamiltonian_matrix(const Atom::AtomicSystem &system, const SplineBasis &basis, int l_number) -> SquareMatrix<double>;

auto solve_hydrogen_like(Atom::AtomicSystem &atom, SplineBasis &basis, int l_number, int max_n_number) -> std::vector<Electron>;


#endif