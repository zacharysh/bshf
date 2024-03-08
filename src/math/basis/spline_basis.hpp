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

#include "../../IO/IO.hpp"


class SplineBasis
{
    // better name?
    //typedef std::vector<std::vector<double> > vector_set;
    
    private:
    auto construct_matrix() -> void;
    auto construct_spline_vectors() -> void;

    public:
    int num_spl;
    int k_spline;
    
    LinearGrid r_grid;

    std::vector<std::vector<double> > bspl;
    std::vector<std::vector<double> > bspl_derivative;

    SquareMatrix<double> Bmatrix;

    SplineBasis(): num_spl(), k_spline(), r_grid(), bspl(), bspl_derivative(), Bmatrix() {};
    SplineBasis(const LinearGrid &grid_, int k_spline, int n_spline);

    auto grid_size() const -> std::size_t { return r_grid.grid_size; }
};


#endif