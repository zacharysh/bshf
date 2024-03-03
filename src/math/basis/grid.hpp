#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>

#include "../../IO/IO.hpp"

class Grid
{
    public:
    const double r0;
    const double r_max;
    std::vector<double> range;
    std::size_t grid_size;

    Grid(double r0_, double r_max_, std::vector<double> range_)
    : r0(r0_), r_max(r_max_), range(range_), grid_size(range_.size()) {};

    auto at(std::size_t i) const -> double { return range.at(i); }
};

class LinearGrid : public Grid
{
    public:
    const double dr;

    LinearGrid(double r0_, double r_max_, std::size_t grid_size_)
    : dr((r_max_ - r0_) / (static_cast<int>(grid_size_) - 1)),
    Grid(r0_, r_max_, std::vector<double>(grid_size_))
    {
        IO::msg::construct<double>("linear grid", {{"r_min", r0_}, {"r_max", r_max_}, {"num_points", grid_size_}});
        
        for(int i = 0; i < static_cast<int>(grid_size_); ++i)
            range.at(i) = i * dr + r0_;

        IO::msg::done();
    }
};

#endif