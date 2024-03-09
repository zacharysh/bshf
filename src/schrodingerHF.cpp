#include <iostream>

#include <vector>

#include <cmath> // sqrt

#include <fstream>

#include <string> // std::string, std::stoi
#include <sstream> // std::string, std::stoi

#include "math/matrix.hpp"
#include "math/math_other.hpp"

#include "math/basis/spline_basis.hpp"
#include "math/basis/grid.hpp"

#include "io.hpp"

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
    
    int Z = 1;
    Potential::Type nuclear_potential = Potential::Type::Unknown;
    Potential::Type interaction_type = Potential::Type::Unknown;


    // std::size_t ?
    int N_grid_size = 0;

    bool fill_atom = true;
    bool gen_spectrum = false;
    bool calc_lifetime = false;

    std::pair<int, int> excited_valence_data {};
    bool excited_valence = false;

    for (int i = 1; i < argc; ++i)
    {
        // Next input parameter.
        std::string param = argv[i];
        std::string next_arg;
        if (argv[i+1])
            next_arg = argv[i+1];
        else break;
        

        if (param == "-Z")
        {
            Z = std::stoi(next_arg);
        }
        else if ((param == "--excited-valence"))
        {
            if(next_arg == "2p") {excited_valence_data = {2, 1}; excited_valence = true; }
        }
        else if ((param == "--nuclear-potential" || param == "-np") && argc > i)
        {
            // tolower?
            if(next_arg == "coulomb" || next_arg == "Coulomb")
                nuclear_potential = Potential::Type::Coulomb;
            else
            {
                IO::log(LogType::error, "Invalid nuclear potential given.", next_arg);
                return EXIT_FAILURE;
            }
        }
        else if ((param == "--interaction-potential"|| param == "-ip"))
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
                IO::log(LogType::error, "Invalid interaction potential given.", next_arg);
                return EXIT_FAILURE;
            }
        }

        else if ((param == "--grid-size"))
            N_grid_size = std::stoi(next_arg);
        
        else if ((param == "--calc-lifetime" || param == "--calculate-lifetime"))
            calc_lifetime = true;
    }

    if(nuclear_potential == Potential::Type::Unknown)
    {
        IO::log(LogType::warn,"No potential given, defaulting to Coulomb");
        nuclear_potential = Potential::Type::Coulomb;
    }
    if(interaction_type == Potential::Type::Unknown)
    {
        IO::log(LogType::warn,"No interaction potential specified. Ignoring");
    }
    if(N_grid_size == 0)
    {
        N_grid_size = 10001;
        IO::log_params(LogType::warn, "No grid size specified. Using default", {{"n", N_grid_size}});
    }
    
    constexpr auto r0 = 5.0e-6;
    constexpr auto rmax = 20.0;
    constexpr auto k_spline = 7;
    constexpr auto n_spline = 90;
    

    SplineBasis basis {LinearGrid(r0, rmax, N_grid_size), k_spline, n_spline};

    Atom atomic_system(Z, nuclear_potential, interaction_type, basis);
    
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
            greens_perturbation(atomic_system, atomic_system.electrons.back());
            IO::done(-1);

            if(excited_valence)
            {
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);

                greens_perturbation(atomic_system, atomic_system.electrons.back());
            }
            break;
        }
        case Potential::Type::HF_Direct:
        {
            HartreeFock::solve_self_consistent(atomic_system);
            IO::done(-1);
            
            if(excited_valence)
            {
                solve_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
            }
            break;
        }
        case Potential::Type::HartreeFock:
        {
            HartreeFock::solve(atomic_system);
            IO::done(-1);

            if(excited_valence)
            {
                HartreeFock::solve_full_excited_valence(atomic_system, excited_valence_data.first, excited_valence_data.second);
            }
                
            break;
        }
        default:
            return EXIT_FAILURE;
    }

    atomic_system.print_states();

    const double expt_2p_energy = -0.13023;
    const double expt_2p_lifetime = 27.102; // ns

    // 2p1/2 lifetime: 1.80415% error
    if(calc_lifetime && excited_valence)
    {
        auto predicted_2p_lifetime = calculate_lifetime(atomic_system.electrons.at(1), atomic_system.electrons.at(2), atomic_system.basis.r_grid);
        predicted_2p_lifetime *= 10e8;
        
        std::cout   << "\n2p energy: " << atomic_system.valence().energy << "au (expt. " << expt_2p_energy << "au, "
                    << abs(abs(atomic_system.electrons.back().energy) - abs(expt_2p_energy)) / abs(expt_2p_energy) * 100.0
                    << "% error).";

        std::cout   << "\n2p lifetime: " << predicted_2p_lifetime << "ns (expt. " << expt_2p_lifetime << "ns, "
                    << abs(abs(predicted_2p_lifetime) - abs(expt_2p_lifetime)) / abs(expt_2p_lifetime) * 100.0
                    << "% error).\n";
    }
    std::ofstream ofs;
    ofs.open("output/Li.txt", std::ofstream::out | std::ofstream::trunc);

    ofs << "r";
    for (int i = 0; i < (int)atomic_system.electrons.size(); ++i)
    {
        ofs << ", " << atomic_system.electrons.at(i).state_label;
    }
    ofs << "\n";
    for(int i = 0; i < N_grid_size; ++i)
    {
        ofs << basis.r_grid.at(i);
        for(int j = 0; j < (int)atomic_system.electrons.size(); ++j)
        {
            ofs << ", " << atomic_system.electrons.at(j).amplitude.at(i);
        }
        ofs << "\n";
    }
    ofs.close();
    

    return EXIT_SUCCESS;
}