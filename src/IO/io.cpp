#include "IO.hpp"

inline
auto IO::msg::print_msg(std::string text) -> void
{
    if(verbose)
        std::cout << text;
}

template <typename T>
inline
auto IO::msg::print_msg(T param) -> void
{
    if(verbose)
        std::cout << param;
}


auto IO::msg::print_msg(std::string action, std::string text) -> void
{
    IO::msg::print_msg(action + " " + text);
}

auto IO::msg::print_new_msg() -> void
{
    for (int i = 0; i < depth; ++i)
        IO::msg::print_msg("  ");
    IO::msg::print_msg("> ");
}

auto IO::msg::print_new_msg(std::string action, std::string name) -> void
{
    print_new_msg();
    print_msg(action, name);
}

auto IO::msg::done(bool up_layer) -> void
{ 
    if(up_layer)
    {
        print_new_msg();
        IO::msg::depth -= 1;
        print_msg("\033[0;32mDone\033[0m.\n");
    }
    else
        print_msg(" \033[0;32mdone\033[0m.\n");
}

auto IO::msg::action(std::string action, std::string name, bool down_layer) -> void
{
    print_new_msg(action, name);
    
    if(down_layer)
    {
        depth += 1;
        print_msg(":\n");
    }
    else
        print_msg("...");
}

template <typename... T>
auto IO::msg::print_values(std::initializer_list<std::pair<std::string, T...> > params) -> void
{
    print_msg(" (");
    
    for (std::pair<std::string, T...> iter : params)
    {
        print_msg("\033[0;36m");
        
        if(iter.first != "")
        {
            print_msg(iter.first, "= ");
        }
        print_msg<T...>(iter.second);
        print_msg("\033[0;0m");

        if(iter != *(params.end()-1))
            print_msg(", ");
    }
    print_msg(")");
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
        IO::msg::print_msg(":\n");
    }
    else
        IO::msg::print_msg("...");
}
// FIX ME!
auto IO::msg::error(std::pair<std::string, int> params, bool down_layer) -> void
{
    IO::msg::print_msg<int>("\033[0;31mError\033[0;0m", {params}, down_layer);
    std::cout << ".\n";
}

auto IO::msg::error_msg(std::string msg) -> void
{
    IO::msg::print_msg("\033[0;31mError\033[0;0m", msg);
    IO::msg::print_msg(".\n");
}

auto IO::msg::warning(std::pair<std::string, int> params, bool down_layer) -> void
{
    IO::msg::print_new_msg();
    IO::msg::print_msg<int>("\033[0;33mWarning\033[0;0m", {params}, down_layer);
    IO::msg::print_msg(".\n");
}

auto IO::msg::warning_msg(std::string msg) -> void
{
    IO::msg::print_new_msg("\033[0;33mWarning\033[0;0m:", msg);
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
template void IO::msg::construct(std::string msg, std::initializer_list<std::pair<std::string, std::string> > params, bool new_layer); 