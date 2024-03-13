#include "atom.hpp"

// Is there a better way to do this?
auto print_state_label(const Electron &psi, bool excited_state_present) -> void
{
    auto label = psi.state_label;

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

auto Atom::print_states() const -> void
{
    // Perhaps we wish to include the normalisation?
    
    const std::string integral("\u222Bdr|P(r)|\u00B2");
    const std::string epsilon("\u03B5");
    const std::string delta("\u0394");

    // Check if we will have to eventually print columns for energy correction.
    bool print_perturbation = false;
    for (auto iter : electrons)
        if(iter.has_correction == true) { print_perturbation = true; break; };
    
    // Construct title.
    std::cout << "\n\033[0;36mResults:\033[0;0m\n\n";
    
    // Construct headings.
    printf(" %7s %14s %14s", "State", (epsilon + " (au)").c_str(), (epsilon + " (eV)").c_str());

    if(print_perturbation == true)
    {
        printf(" %15s %15s", (delta + epsilon + " (au)").c_str(), (delta + epsilon + " (eV)").c_str());
        printf(" %16s %16s", (epsilon + "+" + delta + epsilon + " (au)").c_str(), (epsilon + "+" + delta + epsilon + " (eV)").c_str());
    }

    printf(" %13s  %13s %17s", "<r> (a0)", "<r\u00B2> (a0)", integral.c_str());
    std::cout << "\n";

    bool excited_state_present = false;

    for (const auto &psi : electrons)
    {
        if(psi.filled == false)
            excited_state_present = true;
        // First, print the state label.
        print_state_label(psi, excited_state_present);
        
        // Print energies
        printf(" %13.5f %13.5f", psi.energy, PhysicalConstants::Eh_to_eV(psi.energy));

        if(print_perturbation == true)
        {
            printf(" %13.5f %13.5f", psi.energy_correction, PhysicalConstants::Eh_to_eV(psi.energy_correction));
            printf(" %13.5f %13.5f", psi.energy + psi.energy_correction, PhysicalConstants::Eh_to_eV(psi.energy + psi.energy_correction));
        }

        // Print radial moments
        printf(" %13.5f %13.5f %13.5f\n", psi.calculate_radial_moment(basis.r_grid, 1), psi.calculate_radial_moment(basis.r_grid, 2), simpson_linear(basis.r_grid.dr, psi.P * psi.P));
    }
    std::cout   << (excited_state_present ? "\nKey: '\033[0;96m*\033[0;0m' indicates excited state; '\033[0;33m-\033[0;0m' indicates unfilled state."  : "")
                << (print_perturbation ? ("\nN.B. " + epsilon + " is the unperturbed energy V = V_c + V_Gr. " + delta + epsilon + " is the first-order perturbative correction.") : "")
                << "\n";
}