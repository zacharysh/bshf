#ifndef SPLINE_BASIS_HPP_
#define SPLINE_BASIS_HPP_

#include <vector>
#include <limits>    // std::numeric_limits<double>::epsilon()
#include <algorithm> // std::transform
#include <numeric>   // std::iota
#include <execution> // std::execution

#include "bspline.hpp"
#include "grid.hpp"

#include "../matrix.hpp"

#include "../../io.hpp"


class SplineBasis
{
    // better name?
    //typedef std::vector<std::vector<double> > vector_set;
    
    private:
    auto construct_spline_vectors(double h) -> void;
    auto construct_matrix() -> void;

    public:
    int num_spl;
    int k_spline;
    
    LinearGrid r_grid;

    std::vector<std::vector<double> > bspl;
    std::vector<std::vector<double> > bspl_derivative;

    SquareMatrix<double> Bmatrix;

    SplineBasis(): num_spl(), k_spline(), r_grid(), bspl(), bspl_derivative(), Bmatrix() {};
    SplineBasis(const LinearGrid &grid_, int k_spline_, int n_spline_, double machine_eps_);

    auto grid_size() const -> std::size_t { return r_grid.grid_size; }
};


#endif