#include "io.hpp"

auto IO::print_to_file(std::string file_name, std::vector<std::pair<std::string, std::vector<double>>> &values) -> void
{
    if(print_results == false)
        return;
        
    std::ofstream ofs;

    ofs.open("./output/" + file_name + ".txt", std::ofstream::out | std::ofstream::trunc);


    for (auto iter : values)
    {
        ofs << iter.first;

        // Don't want to print too many commas.
        if(iter != *(values.end() - 1))
            ofs << ", ";
    }
    ofs << "\n";

    for(std::size_t i = 0; i < (values.begin())->second.size(); ++i)
    {
        for (auto iter : values)
        {
            ofs << iter.second.at(i);

            // Don't want to print too many commas.
            if(iter != *(values.end() - 1))
                ofs << ", ";
        }
        ofs << "\n";
    }

    ofs.close();
}