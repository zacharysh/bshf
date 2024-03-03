#include "wavefunction.hpp"

Electron::Electron(const int n_, const int l_, const int m_, const double energy_, const std::vector<double> coeffs_, const SplineBasis &basis_)
: n(n_), l(l_), m(m_),// spin_up(spin_up_),
state_label(std::to_string(n)),
energy(energy_),
coeffs(coeffs_),
basis_size(coeffs.size()),
P(std::vector<double>(basis_.grid_size())),
amplitude(std::vector<double>(basis_.grid_size()))
{
    // Must be a better way
    if( l == 0)
        state_label += "s";
    else if( l == 1)
        state_label += "p";
    else if( l == 2)
        state_label += "d";
    /*
    if(spin_up)
    {
        state += "  " +std::to_string(2*l + 1) + "/2";
    }
    else
    {
        if ((2*l - 1) > 0)
            state += " ";
        state += " " + std::to_string(2*l - 1) + "/2";
    }
    */
        
    //IO::msg::construct<int>("electron state", {{"n", n_}, {"l", l_}, {"m", m_}});
    IO::msg::construct<std::string>("electron state", {{std::string(), state_label}});

    for (std::size_t i = 0; i < basis_.grid_size(); ++i)
    {
        for(std::size_t j = 0; j < basis_size; ++j)
        {
            P.at(i) += basis_.bspl.at(j).at(i) * coeffs.at(j);
        }
        amplitude.at(i) = P.at(i) / basis_.r_grid.range.at(i);
    }

    //IO::msg::print_values<double>({{"E", energy}, {"\u222Bdr|P(r)|\u00B2", trapz_linear(basis.dr, P * P)}});

    r1_moment = calculate_radial_moment(basis_.r_grid, 1);
    r2_moment = calculate_radial_moment(basis_.r_grid, 2);

    IO::msg::print_values<double>({{"\u222Bdr|P(r)|\u00B2", trapz_linear(basis_.r_grid.dr, P * P)}});

    
    IO::msg::done();
}

auto Electron::calculate_radial_moment(const LinearGrid &r_grid, int k) -> double
{   

    auto r_k = r_grid.range;
    // better way?
    for (int i = 1; i < k; ++i)
        r_k = r_k * r_grid.range;

    return trapz_linear(r_grid.dr, P * P * r_k);
}

auto print_states(std::vector<Electron> states) -> void
{
    
    // Construct title.
    // Construct headings.
    
    //\u222Bdr|P(r)|\u00B2

    //std::cout << " n | En    | predicted\n";
    std::cout << "\n\033[0;36mResults:\033[0;0m\n";
    //std::cout << "Key: '\033[0;96m*\033[0;0m' indicates excited state; '\033[0;33m-\033[0;0m' indicates unfilled state.\n";
    std::cout << "Key: '*' indicates excited state; '-' indicates unfilled state.\n";
    printf(" %7s %13s %13s %13s  %13s", "State", "E (au)", "E (eV)", "<r> (a0)", "<r\u00B2> (a0)");
    std::cout << "\n";
    for (auto psi : states)
    {
        printf(" %7s %13.5f %13.5f %13.5f %13.5f\n",
            psi.state_label.c_str(), psi.energy, PhysicalConstants::Eh_to_eV(psi.energy), psi.r1_moment, psi.r2_moment);
    }
}