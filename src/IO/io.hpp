#ifndef IO_HPP_
#define IO_HPP_

#include <string> // std::string

#include <iostream> // std::cout

#include <fstream> // std::fstream

namespace IO
{

    //auto read_file(std::ifstream &input_file) -> std::tuple<double, double, int, int, int>;

namespace msg
{
    auto print_new_msg() -> void;
    auto print_msg(std::string text) -> void;
    auto print_msg(std::string action, std::string text) -> void;
    auto print_new_msg(std::string action, std::string text) -> void;

    inline int depth = 0;

    //is this the best way to do this?


    template <typename... T>
    void print_values(std::initializer_list<std::pair<std::string, T...> > params);

    template <typename... T>
    void print_msg(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

    auto done(bool up_layer = false) -> void;

    template <typename... T>
    void action(std::string action, std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

    auto construct(std::string msg, bool down_layer = false) -> void;
    auto call(std::string msg, bool down_layer = false) -> void;
    void error(std::pair<std::string, int> params, bool down_layer = false);

    template <typename... T>
    void construct(std::string msg, std::initializer_list<std::pair<std::string, T...> > params, bool down_layer = false);

};
}; // namespace IO

#endif