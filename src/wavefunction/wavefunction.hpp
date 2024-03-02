#ifndef WAVEFUNCTION_HPP_
#define WAVEFUNCTION_HPP_

#include <vector>
#include <iostream>

#include "../IO/io.hpp"
#include "../math/SplineBasis.hpp"

#include "../math/math_other.hpp"

class Electron
{
    public:
    int n;
    int l;
    int m;

    double energy;
    std::vector<double> coeffs;
    std::size_t basis_size;

    std::vector<double> P;
    std::vector<double> amplitude;

    //Electron();

    Electron(const int n_, const int l_, const int m_, const double energy_, const std::vector<double> coeffs_, const SplineBasis &basis);
};

#endif