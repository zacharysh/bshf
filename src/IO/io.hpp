#ifndef IO_HPP_
#define IO_HPP_

#include <string> // std::string
#include <iostream> // std::cout
#include <vector> // std::vector

namespace IO
{
    //auto parse_arg_pair(std::initializer_list<std::string> key, std::string arg1, std::string arg2);

namespace msg
{
    auto print_new_msg() -> void;
    auto print_msg(std::string text) -> void;
    auto print_msg(std::string action, std::string text) -> void;
    auto print_new_msg(std::string action, std::string text) -> void;

    template <typename T>
    auto print_msg(T text) -> void;

    inline int depth = 0;
    inline bool verbose = true;

    //is this the best way to do this?

    template <typename... T>
    void print_values(std::initializer_list<std::pair<std::string, T...> > params);

    template <typename... T>
    void print_msg(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

    auto done(bool up_layer = false) -> void;

    template <typename... T>
    void action(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

    auto action(std::string action, std::string msg, bool down_layer = false) -> void;

    auto error_msg(std::string msg) -> void;
    void error(std::pair<std::string, int> params, bool down_layer = false);
    void warning(std::pair<std::string, int> params, bool down_layer = false);
    void warning_msg(std::string msg);

    template <typename... T>
    void construct(std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

};

}; // namespace IO

#endif