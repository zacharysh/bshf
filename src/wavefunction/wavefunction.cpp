#include "wavefunction.hpp"

Electron::Electron(const int n_, const int l_, const int m_, const double energy_, const std::vector<double> coeffs_, const SplineBasis &basis)
: n(n_), l(l_), m(m_),
energy(energy_),
coeffs(coeffs_),
basis_size(coeffs.size()),
P(std::vector<double>(basis.r_grid.size())),
amplitude(std::vector<double>(basis.r_grid.size()))
{
    IO::msg::construct<int>("electron state", {{"n", n_}, {"l", l_}, {"m", m_}});

    for (std::size_t i = 0; i < basis.r_grid.size(); ++i)
    {
        for(std::size_t j = 0; j < basis_size; ++j)
        {
            P.at(i) += basis.bspl.at(j).at(i) * coeffs.at(j);
        }
        amplitude.at(i) = P.at(i) / basis.r_grid.at(i);
    }

    IO::msg::print_values<double>({{"E", energy}, {"\u222Bdr|P(r)|\u00B2 ", trapz_linear(basis.dr, P * P)}});

    
    IO::msg::done();
}