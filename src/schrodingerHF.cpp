#include <iostream>

#include <vector>

#include <cmath> // sqrt

#include <fstream>

#include <string> // std::string, std::stoi
#include <sstream> // std::string, std::stoi

#include <chrono> // std::chrono::high_resolution_clock


#include "io.hpp"

#include "math/matrix.hpp"
#include "math/math_other.hpp"

#include "math/basis/spline_basis.hpp"
#include "math/basis/grid.hpp"

#include "physics/wavefunction.hpp"
#include "physics/atom.hpp"
#include "physics/solve/solve.hpp"
#include "physics/solve/hartree_fock.hpp"
#include "physics/potential.hpp"

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        IO::log(LogType::error, "No input arguments. Aborting");
        return EXIT_FAILURE;
    }

    auto t_start = std::chrono::high_resolution_clock::now();
    
    // Default values for Lithium.
    // TODO: Implement HF Method for all atoms.
    int Z = 3;

    Potential::Type nuclear_potential = Potential::Type::Coulomb;
    Potential::Type interaction_type = Potential::Type::Unknown;

    std::size_t N_grid_size = 2500;

    bool calc_lifetime = false;
    bool excited_valence = false;
    std::pair<int, int> excited_valence_data {};

    int n_max_state = 2;
    int l_max_state = 0;

    bool print_results = false;

    constexpr auto r0 = 1.0e-6;
    constexpr auto rmax = 40.0;
    constexpr auto k_spline = 7;
    constexpr auto n_spline = 60;

    for (int i = 1; i < argc; ++i)
    {
        // Next input parameter.
        std::string param = argv[i];
        std::string next_arg;
        if(param == "--calc-lifetime" || param == "--calculate-lifetime")
        {
            calc_lifetime = true;
            continue;
        }

        if (argv[i+1])
            next_arg = argv[i+1];
        
        if (param == "-Z")
        {
            Z = std::stoi(next_arg);
        }
        else if ((param == "--excited-valence"))
        {
            if(next_arg == "2p")
            {
                excited_valence_data = {2, 1};
                excited_valence = true;
            }
        }
        else if ((param == "--nuclear-potential" || param == "-np") && argc > i)
        {
            if(next_arg == "coulomb" || next_arg == "Coulomb")
                nuclear_potential = Potential::Type::Coulomb;
            else
            {
                IO::log(LogType::error, "Unknown nuclear potential given.", next_arg);
                return EXIT_FAILURE;
            }
        }
        else if (param == "--interaction-potential"|| param == "-ip")
        {
            if(next_arg == "greens" || next_arg == "Greens")
                interaction_type = Potential::Type::Greens;
            else if(next_arg == "self-consistent-hartree" || next_arg == "Self-Consistent-Hartree" || next_arg == "SCH")
            {
                interaction_type = Potential::Type::HF_Direct;
            }
            else if(next_arg == "hartree-fock" || next_arg == "Hartree-Fock" || next_arg == "HF")
            {
                interaction_type = Potential::Type::HartreeFock;
            }
            else
            {
                IO::log(LogType::error, "Unknown interaction potential given.", next_arg);
                return EXIT_FAILURE;
            }
        }
        else if (param == "-v" || param == "--verbose")
        {
            if(next_arg == "false")
                IO::verbose = false;
            else if (next_arg != "true")
            {
                IO::log(LogType::error, "Invalid argument", next_arg);
                return EXIT_FAILURE;
            }
        }
        else if (param == "--print" || param == "--print-to-disc")
        {
            if(next_arg == "false")
                print_results = false;
            else if (next_arg != "true")
            {
                IO::log(LogType::error, "Invalid argument", next_arg);
                return EXIT_FAILURE;
            }
        }
        else if (param == "--grid-size")
            N_grid_size = std::stoi(next_arg);
        else if (param == "--n-max")
            n_max_state = std::stoi(next_arg);
        else if (param == "--l-max")
            l_max_state = std::stoi(next_arg);
    }

    if(nuclear_potential == Potential::Type::Unknown)
    {
        IO::log(LogType::warn,"No potential given, defaulting to Coulomb");
    }
    if(interaction_type == Potential::Type::Unknown)
    {
        IO::log(LogType::warn,"No interaction potential specified. Ignoring");
    }
    // Default size.
    if(N_grid_size == 2501) 
    {
        IO::log_params(LogType::warn, "No grid size specified. Using default", {{"n", N_grid_size}});
    }
    
    SplineBasis basis {LinearGrid(r0, rmax, N_grid_size), k_spline, n_spline};
    Atom atomic_system(Z, n_max_state, l_max_state, nuclear_potential, interaction_type, basis);
    
    // Best way to do this?
    switch (interaction_type)
    {
        case Potential::Type::Unknown:
        {
            solve_atom(atomic_system);
            IO::done(-1);

            if(excited_valence)
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
            break;
        }
        case Potential::Type::Greens:
        {
            solve_atom(atomic_system);
            greens_perturbation(atomic_system, atomic_system.valence());
            IO::done(-1);

            if(excited_valence)
            {
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
                greens_perturbation(atomic_system, atomic_system.valence());
            }
            break;
        }
        case Potential::Type::HF_Direct:
        {
            HartreeFock::solve_self_consistent(atomic_system);
            IO::done(-1);
            
            if(excited_valence)
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
            
            break;
        }
        case Potential::Type::HartreeFock:
        {
            HartreeFock::solve(atomic_system);
            IO::done(-1);

            if(excited_valence)
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
            
            break;
        }
        default:
            return EXIT_FAILURE;
    }

    atomic_system.print_states();


    /************************************************************* //
    //  State    (l_max)     (l_state)    Johnson       Exp.       //
    //  -2s     -0.19631    -0.19631    -0.196304   −0.19814       //
    //  2p      -0.13079    -0.12863    -0.128637   −0.13023       //
    // *************************************************************/

    // Andonis
    // 2s -0.19899
    // 2p -0.128294

    std::cout << "\n";

    // 2p1/2 lifetime: 1.80415% error
    if(calc_lifetime && excited_valence)
    {
        constexpr auto expt_2p_energy = -0.13023; // au
        constexpr auto expt_2p_lifetime = 27.102; // ns

        auto predicted_2p_lifetime = calculate_lifetime(atomic_system.electrons.at(1), atomic_system.valence(), atomic_system.basis.r_grid);
        predicted_2p_lifetime *= 1e9;
        
        std::cout   << "2p energy: " << atomic_system.valence().energy << "au (expt. " << expt_2p_energy << "au, "
                    << abs(abs(atomic_system.valence().energy) - abs(expt_2p_energy)) / abs(expt_2p_energy) * 100.0
                    << "% error).\n";
        
        if(interaction_type == Potential::Type::Greens)
        {
            auto corrected_energy = atomic_system.valence().energy + atomic_system.valence().energy_correction;

            std::cout   << "2p energy (w/ pert.): " << corrected_energy  << "au (expt. " << expt_2p_energy << "au, "
                        << abs(abs(corrected_energy) - abs(expt_2p_energy)) / abs(expt_2p_energy) * 100.0
                        << "% error).\n";
        }

        std::cout   << "2p lifetime: " << predicted_2p_lifetime << "ns (expt. " << expt_2p_lifetime << "ns, "
                    << abs(abs(predicted_2p_lifetime) - abs(expt_2p_lifetime)) / abs(expt_2p_lifetime) * 100.0
                    << "% error).\n";
    }



    // Save output.
    std::cout << "\nPrinting results to disc.";

    // Temp. vector to store all electrons plus radial grid in a reusable format.
    std::vector<std::pair<std::string, std::vector<double> > > print_vector(1 + atomic_system.electrons.size());

    print_vector.front() = {"r", atomic_system.get_range()};

    std::transform(atomic_system.electrons.begin(), atomic_system.electrons.end(), std::next(print_vector.begin()),
        [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.amplitude}; });

    IO::print_to_file("Z=" + std::to_string(atomic_system.Z) + "_amp", print_vector);

    std::transform(atomic_system.electrons.begin(), atomic_system.electrons.end(), std::next(print_vector.begin()),
        [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.P}; });

    IO::print_to_file("Z=" + std::to_string(atomic_system.Z) + "_P", print_vector);

    std::transform(atomic_system.electrons.begin(), atomic_system.electrons.end(), std::next(print_vector.begin()),
        [] (Electron &el) -> std::pair<std::string, std::vector<double> > { return {el.state_label, el.P * el.P}; });

    IO::print_to_file("Z=" + std::to_string(atomic_system.Z) + "_P_sq", print_vector);
    
    // Finished.
    auto t_end = std::chrono::high_resolution_clock::now();

    std::cout << "\nFinished. Total time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << "ms.\n";

    return EXIT_SUCCESS;
}