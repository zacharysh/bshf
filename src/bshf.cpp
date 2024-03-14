#include <iostream>

#include <vector>

#include <cmath> // std::sqrt, std::cbrt

#include <string> // std::string, std::stoi

#include <chrono> // std::chrono::high_resolution_clock

#include "io.hpp"

#include "math/basis/grid.hpp"
#include "math/basis/spline_basis.hpp"

#include "physics/potential.hpp"
#include "physics/atom.hpp"
#include "physics/wavefunction.hpp"

#include "physics/solve/solve.hpp"
#include "physics/solve/hartree_fock.hpp"

auto print_doc() -> void
{
    
    std::cout 
        << "Non-relativistic B-Spline Hartree-Fock atomic many-body solver.\n"
        << "Last updated: " << __TIMESTAMP__ << ", compilation date: " << __DATE__ << ".\n"
        << "Author: Zachary Stevens-Hough.\n"
        << "Usage: 'bshf [options]'.\n"
        << "Example:\n\t'bshf --Z 3 --n-max 10 --l-max 1' computes the first 10 s-wave and p-wave orbitals for hydrogen-like lithium.\n"
        << "Options:\n";

        printf("%+30s:\n%+30s: print this message.\n", "--help", "-h");
        printf("\n%+30s:\n%+30s: print values to disc (saves to './output/') (default: false).\n", "--print-to-disc", "-p");
        printf("\n%+30s:\n%+30s: display verbose output (default: false).\n", "--verbose", "-v");
        
        printf("\n%+30s: set atomic number (default: Z = 3).\n", "--Z [Z]");
        printf("%+30s: set maximum principle quantum number to solve (default: n_max = 2).\n", "--n-max [n_max]");
        printf("%+30s: set maximum angular quantum number to solve (default: l_max = 1).\n", "--l-max [l_max]");
        printf("\n%+30s:\n%+30s: set nuclear potential to 'pot', where 'pot' \u2208 {'coulomb'} (default: 'coulomb').\n", "--nuclear-potential [pot]", "-np [pot]");
        printf("\n%+30s:\n%+30s: set interaction potential to 'pot', where 'pot' \u2208 {'greens', 'hartree', 'hartree-fock'} (default: '').\n", "--interaction-potential [pot]", "-ip [pot]");
}

int main(int argc, char **argv)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    if (argc == 1)
    {
        print_doc();
        return EXIT_SUCCESS;
    }

    /************************************************************************************************************/
    auto io_verbosity = false;

    Potential::Type nuclear_potential = Potential::Type::Unknown;
    Potential::Type interaction_type = Potential::Type::Unknown;

    // Default values for Lithium.
    // TODO: Implement HF Method for all atoms.
    int Z = 3;  
    int n_max_state = 2;
    int l_max_state = 0;

    // In relation to calculating 2p lifetime.
    bool calc_lifetime = false;
    std::pair<int, int> excited_valence_data {};

    // Grid & basis parameters. Grid size is a variable for testing, but the others have been set to the `optimised' values.
    // Ideally, N_grid_size is odd for the purposes of *slightly* better Simpson integration.
    std::size_t N_grid_size = 2501; 
    constexpr auto r0 = 1.0e-6;
    constexpr auto rmax = 40.0;
    constexpr auto k_spline = 7;
    constexpr auto n_spline = 65;
    /************************************************************************************************************/
    
    // Parse input parameters.
    /************************************************************************************************************/
    for (int i = 1; i < argc; ++i)
    {
        // Next input parameter.
        std::string param = argv[i];
        std::string next_arg = argv[i+1] ? argv[i+1] : "";

        // Parameters for PHYS4070 assignment.
        if(param == "--calculate-lifetime")
        {
            excited_valence_data = {2, 1}; // 2p.
            calc_lifetime = true;
            continue;
        }
        if (param == "-v" || param == "--verbose")
        {
            io_verbosity = true;
            continue;
        }
        if (param == "-p" || param == "--print-to-disc")
        {
            IO::print_results = true;
            continue;
        }
        if (param == "-h" || param == "--help")
        {
            print_doc();
            return EXIT_SUCCESS;
        }

        if(next_arg.empty())
            continue;
        
        if (param == "--Z")
            Z = std::stoi(next_arg);
        else if (param == "--grid-size")
            N_grid_size = std::stoi(next_arg);
        else if (param == "--n-max")
            n_max_state = std::stoi(next_arg);
        else if (param == "--l-max")
            l_max_state = std::stoi(next_arg);
        else if ((param == "--nuclear-potential" || param == "-np") && argc > i)
        {
            if(next_arg == "coulomb")
                nuclear_potential = Potential::Type::Coulomb;
            else
            {
                IO::log(LogType::error, "Unknown nuclear potential given.", "'" + next_arg + "'");
                return EXIT_FAILURE;
            }
        }
        else if (param == "--interaction-potential"|| param == "-ip")
        {
            if(next_arg == "greens")
                interaction_type = Potential::Type::Greens;
            else if(next_arg == "hartree" || next_arg == "SCH")
                interaction_type = Potential::Type::HF_Direct;
            else if(next_arg == "hartree-fock" || next_arg == "HF")
                interaction_type = Potential::Type::HartreeFock;
            else
            {
                IO::log(LogType::error, "Unknown interaction potential given.", "'" + next_arg + "'");
                return EXIT_FAILURE;
            }
        }
    }
    
    if(N_grid_size == 2501) 
        IO::log_params(LogType::warn, "No grid size specified. Using default", {{"n", N_grid_size}});

    if(nuclear_potential == Potential::Type::Unknown)
    {
        IO::log(LogType::warn,"No potential given, defaulting to Coulomb");
        nuclear_potential = Potential::Type::Coulomb;
    }
    if(interaction_type == Potential::Type::Unknown)
        IO::log(LogType::warn,"No interaction potential specified. Solving H-Like system");
    
    IO::verbose = io_verbosity;
    /************************************************************************************************************/

    // When solving Hartree-Fock method, we require higher precision near the origin for our 1s orbitals, so take cbrt so there is less
    // numerical error for r << 1 at the cost of decreased relative precision for r >> 1.
    // For other methods, we want a higher level of precision _everywhere_ (but r << 1 is less important),
    // so we should instead take the sqrt. 
    const auto machine_eps = 
        (interaction_type == Potential::Type::HartreeFock || interaction_type == Potential::Type::HF_Direct)
        ? std::cbrt(std::numeric_limits<double>::epsilon()) : std::sqrt(std::numeric_limits<double>::epsilon());

    
    // Initialise basis and atom and then solve the damned thing.
    /************************************************************************************************************/
    SplineBasis basis {LinearGrid {r0, rmax, N_grid_size}, k_spline, n_spline, machine_eps};
    Atom        atom  {Z, n_max_state, l_max_state, nuclear_potential, interaction_type, basis};
    
    IO::log_params(LogType::info, "Constructing atom", {{"Z", Z}}, 1);

    switch (interaction_type)
    {
        case Potential::Type::Unknown:      
        case Potential::Type::Greens:       { generate_atom(atom);      break; }

        case Potential::Type::HF_Direct:    
        case Potential::Type::HartreeFock:  { HartreeFock::solve(atom); break; }

        // Shouldn't be possible to set the int. potential to V_ee = -Z / r.
        case Potential::Type::Coulomb:      { std::cerr << "How did this even occur!?"; return EXIT_FAILURE; }
    }
    if(calc_lifetime)
        solve_excited_valence(atom, excited_valence_data.first, excited_valence_data.second);
    
    IO::done(-1);
    /************************************************************************************************************/

    atom.print_states();

    // Entirely for the calculations related to the Lithium 2p -> 2s lifetime. See 4070 task sheet.
    /************************************************************************************************************/
    if(calc_lifetime)
    {
        constexpr auto expt_2p_energy = -0.13023; // au
        constexpr auto expt_2p_lifetime = 27.102; // ns

        auto predicted_2p_lifetime = calculate_lifetime(atom.electrons.at(1), atom.valence(), atom.basis.r_grid) * 1e9;
        
        std::cout   << "\n2p energy: " << atom.valence().energy << "au (expt. " << expt_2p_energy << "au, "
                    << abs(abs(atom.valence().energy) - abs(expt_2p_energy)) / abs(expt_2p_energy) * 100.0
                    << "% error).";
        
        if(interaction_type == Potential::Type::Greens)
        {
            auto corrected_energy = atom.valence().energy + atom.valence().energy_correction;

            std::cout   << "\n2p energy (w/ pert.): " << corrected_energy  << "au (expt. " << expt_2p_energy << "au, "
                        << abs(abs(corrected_energy) - abs(expt_2p_energy)) / abs(expt_2p_energy) * 100.0
                        << "% error).\n";
        }

        std::cout   << "\n2p lifetime: " << predicted_2p_lifetime << "ns (expt. " << expt_2p_lifetime << "ns, "
                    << abs(abs(predicted_2p_lifetime) - abs(expt_2p_lifetime)) / abs(expt_2p_lifetime) * 100.0
                    << "% error).\n";
    }
    /************************************************************************************************************/


    // Save output.
    /************************************************************************************************************/
    if(IO::print_results == true)
    {
        std::cout << "\nPrinting results to disc.\n";

        // Temp. vector to store all electrons plus radial grid in a reusable format.
        std::vector<std::pair<std::string, std::vector<double> > > print_vector(1 + atom.electrons.size());

        print_vector.front() = {"r", atom.get_range()};

        // Print amplitudes.
        std::transform(atom.electrons.begin(), atom.electrons.end(), std::next(print_vector.begin()),
            [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.amplitude}; });
        IO::print_to_file("Z=" + std::to_string(atom.Z) + "_amp", print_vector);

        // Print P.
        std::transform(atom.electrons.begin(), atom.electrons.end(), std::next(print_vector.begin()),
            [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.P}; });
        IO::print_to_file("Z=" + std::to_string(atom.Z) + "_P", print_vector);


        // Print |P|^2.
        std::transform(atom.electrons.begin(), atom.electrons.end(), std::next(print_vector.begin()),
            [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.P * el.P}; });
        IO::print_to_file("Z=" + std::to_string(atom.Z) + "_P_sq", print_vector);
    }
    /************************************************************************************************************/


    // Finished.
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout   << "\nFinished. Total time elapsed = "
                << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms.\n";

    return EXIT_SUCCESS;
}