#ifndef SPLINE_BASIS_HPP_
#define SPLINE_BASIS_HPP_

#include "bspline.hpp"

#include "../IO/io.hpp"

#include "matrix.hpp"

#include <vector>

class SplineBasis
{
    typedef std::vector<std::vector<double> > vector_set;
    
    private:
    auto construct_matrix() -> void;
    auto construct_spline_vectors() -> void;

    public:
    int num_spl;
    int k_spline;
    std::vector<double> r_grid;
    double dr;
    vector_set bspl;
    vector_set bspl_derivative;

    SquareMatrix<double> Bmatrix;

    SplineBasis(const std::vector<double> &r_grid_, int k_spline, int n_spline);
};


#endif