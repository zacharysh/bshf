#include <iostream>

#include <vector>

#include <cmath> // sqrt

#include <fstream>

#include <string> // std::string, std::stoi
#include <sstream> // std::string, std::stoi

#include "IO/io.hpp"

#include "wavefunction/wavefunction.hpp"
#include "wavefunction/atom.hpp"
#include "wavefunction/solve.hpp"

#include "math/matrix.hpp"
#include "math/math_other.hpp"

/*
if (param == "--input-file")
{
    // check that theres at least enough arguments for the file name.
    if(argc > i) 
    {
        std::ifstream input_file;
        std::string if_name {argv[i+1]};
        input_file.open(if_name);
        // check successful reading
        if(input_file.is_open())
        {
            const auto [r0, rmax, k_spline, n_spline, N_grid] = IO::read_file(input_file);
        }
    }
    else
    {
        std::cerr << "Not enough input arguments. Needs file name.\n";
        return EXIT_FAILURE;
    }

    
    auto pos = param.find("-Z");
    
    // probably wont work for multiple input arguments.
    if(pos != std::string::npos)
    {
        Z = std::stoi(param.substr(pos + 3));
    }
    else
    {
        std::cout << "WARNING: defaulting to Z = 1.\n";
    }
    
*/

int main(int argc, char **argv)
{
    
    if (argc == 0)
    {
        std::cerr << "Incorrect number of args.\n Requires \'-Z\'.'\n";
        return EXIT_FAILURE;
    }
    
    int Z = -1;
    std::vector<int> l{};

    for (int i = 1; i < argc; ++i)
    {
        // Next input parameter.
        std::string param{argv[i]};

        if (param == "-Z" && argc > i)
        {
            Z = std::stoi(argv[i+1]);
        }
        else if ((param == "-l" || param == "-L") && argc > i)
        {
            std::istringstream next_arg_ss(argv[i+1]);

            std::string token;
            while(getline(next_arg_ss, token, ' '))
            {
                if(std::isdigit(token[0]))
                    l.push_back(std::stoi(token));
            }
        }
    }


    const double r0 = 1.0e-10;
    const double rmax = 25.0;
    const int k_spline = 7;
    const int n_spline = 60;
    const int N_grid = 20000;
    
    std::vector<double> r_grid = construct_grid_linear(r0, rmax, N_grid);
    

    SplineBasis basis(r_grid, k_spline, n_spline);
    Atom::AtomicSystem Li(Z, r_grid);
    
    std::vector<Electron> electrons{};
    for (auto iter : l)
    {
        auto sols = solve_hydrogen_like(Li, basis, iter, 5);
        electrons.reserve(electrons.size() + sols.size());
        electrons.insert(electrons.end(), sols.begin(), sols.end());
    }

    std::vector<double> predicted(electrons.size());
    for (int i = 1; i < (int)electrons.size(); ++i) 
        predicted.at(i-1) = -Z * Z / (2.0 * i * i);

    /*
    std::cout << "---- Energies ----" << '\n';
    std::cout << " n | En    | predicted\n"; 
    for (int i=1; i <= k; ++i)
    {
        printf("%2d | %1.4f | %1.4f\n", i, energies.at(i-1), predicted.at(i-1));
    }
    */
    
    std::ofstream ofs;
    ofs.open("output/Li.txt", std::ofstream::out | std::ofstream::app);

    ofs << "r";
    for (int i = 0; i < (int)electrons.size(); ++i)
    {
        ofs << ", " << i + 1;
    }
    ofs << "\n";
    for(int i = 0; i < N_grid; ++i)
    {
        ofs << r_grid.at(i);
        for(int j = 0; j < (int)electrons.size(); ++j)
        {
            ofs << ", " << electrons.at(j).P.at(i);
        }
        ofs << "\n";
    }
    ofs.close();

    
    return EXIT_SUCCESS;
}