#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>

#include <numeric> // std::iota
#include <algorithm> // std::transform

#include "../../io.hpp"

class LinearGrid
{
    public:
    std::size_t grid_size;
    double r0;
    double r_max;
    double dr;
    std::vector<double> range;
    std::vector<double> range_inv;


    LinearGrid() : r0(), r_max(), range(), grid_size() {};

    LinearGrid(double r0_, double r_max_, std::size_t grid_size_)
    : grid_size(grid_size_),
    r0(r0_), r_max(r_max_),
    dr((r_max_ - r0_) / (static_cast<int>(grid_size_) - 1.0)),
    range(std::vector<double>(grid_size_)),
    range_inv(std::vector<double>(grid_size_))
    {
        IO::log_params(LogType::info, "Constructing linear grid", {{"r_min", r0_}, {"r_max", r_max_}, {"num_points", grid_size_}});

        std::generate(range.begin(), range.end(),
            [n=0, *this] () mutable
                { return std::fma(n++,  dr, r0); });
        

        std::transform(range.begin(), range.end(), range_inv.begin(),
            [] (auto r) { return 1.0 / r; });

        IO::done();
    };

    inline auto at(std::size_t i) const -> double { return range.at(i); }
    inline auto inv_at(std::size_t i) const -> double { return range_inv.at(i); }
};

#endif