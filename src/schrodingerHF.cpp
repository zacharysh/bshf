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

#include "IO/io.hpp"

#include "physics/wavefunction.hpp"
#include "physics/atom.hpp"
#include "physics/solve/solve.hpp"
#include "physics/potential.hpp"

int main(int argc, char **argv)
{
    
    if (argc == 1)
    {
        IO::msg::error_msg("No input arguments. Aborting");
        return EXIT_FAILURE;
    }
    
    int Z = 1;
    Potential::Type nuclear_potential = Potential::Type::Unknown;
    Potential::Type interaction_type = Potential::Type::Unknown;

    std::vector<int> l_vector {};

    // std::size_t ?
    int N_grid_size = 0;

    bool fill_atom = true;
    bool gen_spectrum = false;

    std::pair<int, int> excited_valence_data {};
    bool excited_valence = false;

    for (int i = 1; i < argc; ++i)
    {
        // Next input parameter.
        std::string param{argv[i]};
        

        if (param == "-Z" && argc > i)
        {
            Z = std::stoi(argv[i+1]);
        }
        if ((param == "-L" || param == "-l") && argc > i)
        {
            std::stringstream next_arg_ss(argv[i+1]);

            std::string token;
            while(getline(next_arg_ss, token, ' '))
            {
                if(std::isdigit(token[0]))
                    l_vector.push_back(std::stoi(token));
            }
        }

        if ((param == "--excited-valence") && argc > i)
        {
            // Fix me.
            std::string state_arg(argv[i+1]);
            if(state_arg == "2p")
                {excited_valence_data = {2, 1}; excited_valence = true; }
        }
                
        // Get rid of this?
        else if ((param == "--nuclear-potential" || param == "-np") && argc > i)
        {
            std::string next_arg(argv[i+1]);

            // tolower?
            if(next_arg == "coulomb" || next_arg == "Coulomb")
                nuclear_potential = Potential::Type::Coulomb;
            else
            {
                std::cout << "Invalid nuclear potential '" << next_arg << "' given.\n";
                return EXIT_FAILURE;
            }
        }

        else if ((param == "--interaction-potential"|| param == "-ip") && argc > i)
        {
            std::string next_arg(argv[i+1]);

            // tolower?
            if(next_arg == "greens" || next_arg == "Greens")
                interaction_type = Potential::Type::Greens;
            else if(next_arg == "hartree-fock" || next_arg == "Hartree-Fock" || next_arg == "HF")
                interaction_type = Potential::Type::HartreeFock;
            else
            {
                std::cout << "Invalid nuclear potential '" << next_arg << "' given.\n";
                return EXIT_FAILURE;
            }
        }

        else if ((param == "--grid-size") && argc > i)
        {
            //std::istringstream next_arg_ss(argv[i+1]);
            //next_arg_ss >> N_grid_points;
            N_grid_size = std::stoi(argv[i+1]);
        }
        else if ((param == "--fill-atom"))
        {
            fill_atom = true;
        }
        else if ((param == "--generate-spectrum"))
        {
            gen_spectrum = true;
            fill_atom = false;
        }
    }

    if(nuclear_potential == Potential::Type::Unknown)
    {
        IO::msg::warning("No potential given, defaulting to Coulomb");
        nuclear_potential = Potential::Type::Coulomb;
    }
    if(interaction_type == Potential::Type::Unknown)
    {
        IO::msg::warning("No interaction potential specified. Ignoring");
    }
    if(N_grid_size == 0)
    {
        IO::msg::warning("No grid size specified. Using default");
        N_grid_size = 10000;
    }
    if(l_vector.empty() && gen_spectrum == true)
    {
        IO::msg::error_msg("No l number given");
        return EXIT_FAILURE;
    }

    const double r0 = 1.0e-8;
    const double rmax = 40.0;
    const int k_spline = 7;
    const int n_spline = 60;
    

    SplineBasis basis {LinearGrid(r0, rmax, N_grid_size), k_spline, n_spline};

    Atom atomic_system(Z, nuclear_potential, interaction_type, basis);
    
    if(interaction_type != Potential::Type::HartreeFock)
    {
        // Best way to do this?

        if(fill_atom == true)
            solve_atom(atomic_system);
        else if(gen_spectrum == true)
            for (auto iter : l_vector)
                solve_schrodinger(atomic_system, iter, false);
        
        if(excited_valence)
        {
            // In an atom such as lithium, the valence electron is the only state occupying the valence shell.
            // e.g. lithium first excited state is 1s2 2s1 -> 1s2 2p1).

            //atomic_system.electrons.pop_back();
            //atomic_system.electrons.back().state_label = "\033[0;33m-" + atomic_system.electrons.back().state_label + "\033[0;0m";
            atomic_system.electrons.back().state_label.insert(0, "-");
            atomic_system.electrons.back().filled = false;

            solve_schrodinger(atomic_system, excited_valence_data.second, excited_valence_data.first);

            //atomic_system.electrons.back().state_label = "\033[0;96m*" + atomic_system.electrons.back().state_label + "\033[0;0m";
            atomic_system.electrons.back().state_label.insert(0, "*");
        }
    }
    else
    {

    }

    print_states(atomic_system.electrons);
    
    
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
            ofs << ", " << atomic_system.electrons.at(j).P.at(i);
        }
        ofs << "\n";
    }
    ofs.close();
    

    return EXIT_SUCCESS;
}