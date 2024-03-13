#include "spline_basis.hpp"


SplineBasis::SplineBasis(const LinearGrid &r_grid_, int k_spline_, int n_spline_)
:   num_spl(n_spline_ - 3),
    k_spline(k_spline_),
    r_grid(r_grid_),
    bspl(std::vector<std::vector<double> >(num_spl, std::vector<double>(r_grid.grid_size))),
    bspl_derivative(std::vector<std::vector<double> >(num_spl, std::vector<double>(r_grid.grid_size))),
    Bmatrix(SquareMatrix<double>(num_spl))
{
    IO::log_params(LogType::info, "Constructing B-Spline basis", {{"k", k_spline}, {"n", num_spl}}, 1);

    construct_spline_vectors();
    construct_matrix();

    IO::done(-1);
}

auto SplineBasis::construct_spline_vectors() -> void
{
    IO::log("Constructing basis vectors");
    
    BSpline splines(k_spline, num_spl + 3, r_grid.r0, r_grid.r_max);

    const auto h = std::cbrt(std::numeric_limits<double>::epsilon());
    //const int radix =  std::numeric_limits<double>::radix;
    //const int mantissa = std::numeric_limits<double>::digits;
    //const double h = pow(1.0/radix, mantissa/3);


    for (auto n = 2; n < num_spl + 2; ++n)
    {
        std::transform(r_grid.range.begin(), r_grid.range.end(), bspl.at(n - 2).begin(),
                [splines, n] (auto r) { return splines.b(n, r); });

        std::transform(r_grid.range.begin(), r_grid.range.end(), bspl_derivative.at(n - 2).begin(),
                [h, splines, n] (auto r) { return (splines.b(n, r + h) - splines.b(n, r - h)) / (2*h); });
                //[h, splines, n] (auto r) { return (splines.b(n, r + h) - splines.b(n, r - h)) / ((r+h)-(r-h)); });
                //[h, splines, n] (auto r) { return (splines.b(n, r + h * r) - splines.b(n, r - h * r)) / ((r + h * r) - (r - h * r)); });
    }

    IO::done();
}

auto SplineBasis::construct_matrix() -> void
{
    IO::log("Constructing B-matrix");

    for(int i = 0; i < num_spl; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            Bmatrix(i,j) = simpson_linear(r_grid.dr, bspl.at(i) * bspl.at(j));
        }
    }

    IO::done();
}