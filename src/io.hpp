#ifndef IO_HPP_
#define IO_HPP_

#include <string> // std::string
#include <sstream> // std::string
#include <iostream> // std::cout

#include <vector>
#include <fstream>

#define RED     "\033[0;31m"
#define GREEN   "\033[0;32m"
#define YELLOW  "\033[0;33m"
#define BLUE    "\033[0;36m"
#define CLEAR   "\033[0;0m"

enum class LogType
{
    text,
    info,
    warn,
    error,
    done
};

namespace IO
{

inline int  depth    = 0;
inline bool verbose  = true;

inline bool print_results = false;

inline
auto log(const std::string_view msg, int delta_depth = 0) -> void
{
    if(!verbose)
        return;

    for (int i = 0; i < depth; ++i)
        std::cout << "  ";

    std::cout << "> " << msg;

    depth += delta_depth;

    
    if(delta_depth > 0)
        std::cout << ":\n";
    else
        std::cout << "... ";
}

inline
auto log(LogType msg_type, const std::string_view msg, int delta_depth = 0) -> void
{
    if(!verbose)
        return;

    if(msg_type != LogType::done || (msg_type == LogType::done && delta_depth < 0))
    {
    for (int i = 0; i < depth; ++i)
        std::cout << "  ";

    std::cout << "> ";
    }

    depth += delta_depth;

    switch(msg_type)
    {
        case LogType::text:  { std::cout                           << msg;             break;};
        case LogType::info:  { std::cout  << BLUE                  << msg << CLEAR;    break;};
        case LogType::warn:  { std::cout  << YELLOW << "Warning: " << msg << CLEAR;    break;};
        case LogType::error: { std::cout  << RED    << "Error: "   << msg << CLEAR;    break;};
        case LogType::done:  { std::cout  << GREEN                 << msg << CLEAR;    break;};
    }
    
    if(msg_type != LogType::done && msg_type != LogType::warn && msg_type != LogType::error)
    {
        if(delta_depth > 0)
            std::cout << ":\n";
        else
            std::cout << "... ";
    }
    else
        std::cout << ".\n";
}

inline
auto log(LogType msg_type, const std::string_view msg, const std::string_view info, int delta_depth = 0) -> void
{
    if(!verbose)
        return;


    if(msg_type != LogType::done)
    {
    for (int i = 0; i < depth; ++i)
        std::cout << "  ";

    std::cout << "> ";
    }

    depth += delta_depth;

    switch(msg_type)
    {
        default:             { std::cout                           << msg;             break;};
        case LogType::warn:  { std::cout  << YELLOW << "Warning: " << msg << CLEAR;    break;};
        case LogType::error: { std::cout  << RED    << "Error: "   << msg << CLEAR;    break;};
    }
    switch(msg_type)
    {
        default: break;
        case LogType::info:  { std::cout  << BLUE   << " (" << info << ")" << CLEAR;   break;};
        case LogType::warn:  { std::cout  << YELLOW << " (" << info << ")." << CLEAR;   break;};
        case LogType::error: { std::cout  << RED    << " (" << info << ")." << CLEAR;   break;};
    }

    if(msg_type != LogType::done && msg_type != LogType::warn && msg_type != LogType::error)
    {
    if(delta_depth > 0)
        std::cout << ":\n";
    else
        std::cout << "... ";
    }
    else
        std::cout << "\n";
}


template <typename T = int>
inline
auto log_params(LogType msg_type, const std::string_view msg, const std::initializer_list<T> &params, int delta_depth = 0) -> void
{
    std::stringstream text {""};

    for (const auto &iter : params)
    {
        text << iter;

        if(iter != *(params.end()-1))
            text << ", ";
    }
    
    text << "";

    log(msg_type, msg, text.str(), delta_depth);
}
    
template <typename T = int> 
inline
auto log_params(LogType msg_type, const std::string_view msg, const std::initializer_list<std::pair<std::string, T> > &params, int delta_depth = 0) -> void
{
    std::stringstream text {""};

    for (const std::pair<std::string, T> iter : params)
    {
        text << iter.first << " = " << iter.second;

        if(iter != *(params.end()-1))
            text << ", ";
    }
    
    log(msg_type, msg, text.str(), delta_depth);
}

inline
auto done(const int delta_depth = 0) -> void
{
    if(delta_depth < 0) log(LogType::done, "Done", delta_depth);
    else                log(LogType::done, "done", delta_depth);
}

auto print_to_file(std::string file_name, std::vector<std::pair<std::string, std::vector<double>>> &values) -> void;


}; // namespace IO
#endif