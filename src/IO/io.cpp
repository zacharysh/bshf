#include "io.hpp"

#include <sstream>

auto IO::msg::print_msg(std::string text) -> void
{
    std::cout << text;
}
auto IO::msg::print_msg(std::string action, std::string text) -> void
{
    std::cout << action << " " << text;
}

auto IO::msg::print_new_msg() -> void
{
    for (int i = 0; i < depth; ++i)
        std::cout << "  ";
    std::cout << "> ";
}

auto IO::msg::print_new_msg(std::string action, std::string name) -> void
{
    for (int i = 0; i < depth; ++i)
        std::cout << "  ";
    std::cout << "> ";
    IO::msg::print_msg(action, name);
}

auto IO::msg::done(bool up_layer) -> void
{ 
    if(up_layer)
    {
        print_new_msg();
        IO::msg::depth -= 1;
        std::cout << "\033[0;32mDone\033[0m.\n";
    }
    else
        std::cout << " \033[0;32mdone\033[0m.\n";
}

auto IO::msg::construct(std::string name, bool down_layer) -> void
{
    IO::msg::print_new_msg("Constructing", name);
    
    if(down_layer)
    {
        IO::msg::depth += 1;
        std::cout << ":\n";
    }
    else
        std::cout << "...";
}

auto IO::msg::call(std::string name, bool down_layer) -> void
{
    IO::msg::print_new_msg("Calling", name);

    if(down_layer)
    {
        IO::msg::depth += 1;
        std::cout << ":\n";
    }
    else
        std::cout << "...";
}

template <typename... T>
auto IO::msg::print_values(std::initializer_list<std::pair<std::string, T...> > params) -> void
{
    std::cout << " (";
    
    for (std::pair<std::string, T...> iter : params)
    {
        std::cout << "\033[0;36m" << iter.first << " = " << iter.second << "\033[0;0m";

        if(iter != *(params.end()-1))
            std::cout << ", ";
    }
    std::cout << ")";
}

template <typename... T>
auto IO::msg::print_msg(std::string action, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer) -> void
{
    IO::msg::print_msg(action);
    IO::msg::print_values<T...>(params);
}


template <typename... T>
auto IO::msg::print_msg(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer) -> void
{
    IO::msg::print_msg(action, msg);
    IO::msg::print_values<T...>(params);
}

template <typename... T>
auto IO::msg::action(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer) -> void
{
    IO::msg::print_new_msg();
    IO::msg::print_msg<T...>(action, msg, params);
    
    if(down_layer)
    {
        IO::msg::depth += 1;
        std::cout << ":\n";
    }
    else
        std::cout << "...";
}
// FIX ME!
auto IO::msg::error(std::pair<std::string, int> params, bool down_layer) -> void
{
    IO::msg::print_msg<int>("\033[0;31mError\033[0;0m", {params}, down_layer);
    std::cout << ".\n";
}

template <typename... T>
auto IO::msg::construct(std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer) -> void
{
    IO::msg::action<T...>("Constructing", msg, params, down_layer);
}

template void IO::msg::action(std::string action, std::string msg, std::initializer_list<std::pair<std::string, int> > params, bool new_layer);
template void IO::msg::construct(std::string msg, std::initializer_list<std::pair<std::string, int> > params, bool new_layer); 
template void IO::msg::construct(std::string msg, std::initializer_list<std::pair<std::string, double> > params, bool new_layer); 