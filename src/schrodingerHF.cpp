#include <iostream>

#include <vector>

#include <fstream>

#include "wavefunction/wavefunction.hpp"
#include "wavefunction/atom.hpp"
#include "wavefunction/solve.hpp"

#include "math/matrix.hpp"
#include "math/math_other.hpp"

int main(int argc, char **argv)
{
    /*
    if (argc != 4)
    {
        std::cerr << "Incorrect number of args. Requires 3: N, tmax, v0\n";
        return EXIT_FAILURE;
    }
    */
    
    
    double r_min = 1.0e-5;
    double r_max = 50.0;
    int k_spline = 7;
    int n_spline = 43;

    int N_grid = 10000;

    std::cout << "> Constructing linear grid: r0 = " << r_min << ", r_max = " << r_max << ", N_grid = " << N_grid << "... ";
    std::vector<double> r_grid = construct_grid_linear(r_min, r_max, N_grid);
    std::cout << "done.\n";

    SplineBasis basis(r_grid, k_spline, n_spline);
    Atom::AtomicSystem Li(3, r_grid);
    Wavefunction psi(1, 0, 0); // H-like 1s orbital

    solveHydrogenlike(Li, psi, basis);

    
    std::ofstream ofs;
    ofs.open("1s_test.txt", std::ofstream::out | std::ofstream::app);

    ofs << "r, psi(r)\n";

    for(int i = 0; i < N_grid; ++i)
    {
        ofs << r_grid.at(i) << ", " << psi.amplitude.at(i) << '\n';
        
    }

    ofs.close();
    


    /*
    SquareMatrix<double> sm(2);
    //sm(0,1) = sm(1,0) = 1.0;
    
    for (std::size_t i = 0; i < sm.size_x; ++i) 
    {
        for (std::size_t j = 0; j < sm.size_y; ++j) 
        {
            sm(i, j) = 1.0 / ((int)i + (int)j + 1);
        }
    }
    
    SquareMatrix<double> identity(2);
    for (std::size_t i = 0; i < identity.size_x; ++i) 
    {
        identity(i,i) = 1.0;
    }

  // RealSymmetric: function defined in eigen.hpp
  // MatrixAndVector is a struct defined in eigen.hpp
  // It contians .vector (a std::vector of eigen values)
  // and .matrix (a Matrix of eigenvectors)

  // in c++17, we can do this:
    //const auto [EVectors, EValues] =  
    auto EValues = MatrixTools::solveEigenSystem(sm, identity);
  // otherwise, same as this:
  // eigen::MatrixAndVector temp = eigen::RealSymmetric(sm);
  // const auto &EVectors = tmp.matrix;
  // const auto &EValues  = tmp.vector;


  // Print out the eigenvalues and eigenvectors:
    std::cout << "Eigenvalues: ";
    for (auto v : EValues) {
        std::cout << v << ", ";
    }
    std::cout << '\n';
    std::cout << "With corresponding eigenvectors: \n";
    for (std::size_t i = 0; i < sm.size_x; ++i)
    {
        for (std::size_t j = 0; j < sm.size_y; ++j)
        {
        std::cout << sm(i, j) << ", ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    */
    
    return EXIT_SUCCESS;
}