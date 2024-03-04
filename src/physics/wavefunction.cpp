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
        
    //IO::msg::construct<int>("electron state", {{"n", n_}, {"l", l_}, {"m", m_}});
    IO::msg::construct<std::string>("electron state", {{std::string(), state_label}});

    // Change to scalar multiplication of a vector
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

// We assume a is the only empty state and b is the excited state.
auto calculate_lifetime(Electron a, Electron b, const LinearGrid &r_grid) -> double
{   
    // Use experimental omega.
    const auto omega = 0.06791;

    assert((a.l == 0 && b.l == 1 && a.n == 2 && b.n == 2) && "Only works for 2p -> 2s transition currently.");

    auto Rab = trapz_linear(r_grid.dr, a.P * r_grid.range * b.P);

    auto gamma = 2.0 * Rab * Rab * omega * omega * omega / 3.0 * 1.071e10;

    return 1/gamma;
}

#define EPSILON \u03B5
#define DELTA \u0394

auto print_states(std::vector<Electron> states) -> void
{
    // TODO: Perhaps we wish to include the normalisation?
    //\u222Bdr|P(r)|\u00B2

    std::string epsilon("\u03B5");
    std::string delta("\u0394");

    // Check if we will have to eventually print columns for energy correction.
    bool print_perturbation = false;
    for (auto iter : states)
        if(iter.has_correction == true) { print_perturbation = true; };
    
    // Construct title.
    std::cout << "\n\033[0;36mResults:\033[0;0m\n\n";
    
    // Construct headings.
    printf(" %7s %14s %14s", "State", (epsilon + " (au)").c_str(), (epsilon + " (eV)").c_str());

    if(print_perturbation == true)
    {
        printf(" %15s %15s", (delta + epsilon + " (au)").c_str(), (delta + epsilon + " (eV)").c_str());
        printf(" %16s %16s", (epsilon + "+" + delta + epsilon + " (au)").c_str(), (epsilon + "+" + delta + epsilon + " (eV)").c_str());
    }

    printf(" %13s  %13s", "<r> (a0)", "<r\u00B2> (a0)");
    std::cout << "\n";

    // Bad name...
    bool excited_state_present = false;
    for (const auto &psi : states)
    {
        // First, print the state label.
        print_state_label(psi, excited_state_present);
        
        printf(" %13.5f %13.5f", psi.energy, PhysicalConstants::Eh_to_eV(psi.energy));

        if(print_perturbation == true)
        {
            printf(" %13.5f %13.5f", psi.energy_correction, PhysicalConstants::Eh_to_eV(psi.energy_correction));
            printf(" %13.5f %13.5f", psi.energy + psi.energy_correction, PhysicalConstants::Eh_to_eV(psi.energy + psi.energy_correction));
        }
        printf(" %13.5f %13.5f\n", psi.r1_moment, psi.r2_moment);
    }
    std::cout   << (excited_state_present ? "\nKey: '\033[0;96m*\033[0;0m' indicates excited state; '\033[0;33m-\033[0;0m' indicates unfilled state."  : "")
                << (print_perturbation ? ("\nN.B. " + epsilon + " is the unperturbed energy V = V_c + V_Gr. " + delta + epsilon + " is the first-order perturbative correction.") : "")
                << "\n";
}

// Is there a better way to do this?
inline
auto print_state_label(const Electron &psi, bool &excited_state_present) -> void
{
    auto label = psi.state_label;

    //if(psi.has_correction)
    //    label.push_back('*');
    

    if(psi.filled == false)
        {
            std::cout << "\033[0;33m";
            label = "-" + label;
            printf(" %7s", label.c_str());
            std::cout << "\033[0;0m";

            excited_state_present = true;
        }
        else if(excited_state_present == true)
        {
            std::cout << "\033[0;96m";
            label = "*" + label;
            printf(" %7s", label.c_str());
            std::cout << "\033[0;0m";
        }
        else
            printf(" %7s", label.c_str());
}