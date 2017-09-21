/// Part of: https://github.com/mateidavid/hpptools
/// Commit: 5a6f39c

/// @author    Matei David, Ontario Institute for Cancer Research
/// @version   1.0
/// @date      2013-2017
/// @copyright MIT Public License
///
/// Logger mechanism.
///
/// Properties:
/// - thread-safe, non-garbled output (uses c++11's thread_local)
/// - customizable ostream sink. by default, uses std::clog
///
/// Exports:
/// - macro: LOG (takes 1, 2, or 3 arguments, see below)
/// - macros: LOG_EXIT_, LOG_ABORT, LOG_EXIT, LOG_THROW_, and LOG_THROW
/// - namespace logger
/// - enum logger::level
/// - class logger::Logger
///
/// To use:
/// - In source code, use:
///
///     LOG(info) << "hello" << endl;
///     // or 
///     LOG("main", info) << "hello" << endl;
///     // or 
///     LOG("main", info, sink_os) << "hello" << endl;
///
///   Here, "main" is the facility (a string) and info is the message level.
///   Note that "logger" is a macro which knows how to look up the name info
///   inside logger namespace. The macro introduces C++ code equivalent to:
///
///     if (...message should be ignored...) then; else sink_os
///
///   NOTE: As with assert(), the code in the output stream following the
///   LOG() macro will ***not be executed*** if the log level of the
///   facility is higher than the level of the message.
///
/// - To set the default log level (for unspecified facilities):
///
///     logger::Logger::set_default_level(logger::Logger::level_from_string(s));
///
/// - To set the log level for the "main" facility:
///
///     logger::Logger::set_facility_level("main", logger::Logger::level_from_string(s));
///
/// - To parse a log facility level setting in the form "[<facility>:]<level>":
///
///     logger::Logger::set_level_from_option("alt:debug1", &cerr);
///
/// - By using these functions, one can set log levels using command-line
///   parameters and achieve dynamic log level settings without recompiling.
///
/// - The macros LOG_EXIT_, LOG_ABORT, LOG_EXIT, LOG_THROW_, and LOG_THROW
///   provide a way to specify what to do after logging the message.
///
///     LOG_EXIT_(exit_code)  << "print this message to std::cerr, call std::exit(exit_code)";
///     LOG_ABORT             << "print this message to std::cerr, call std::abort()";
///     LOG_EXIT              << "print this message to std::cerr, call std::exit(EXIT_FAILURE)";
///     LOG_THROW_(Exception) << "throw Exception with this message";
///     LOG_THROW             << "throw std::runtime_error with this message";

#ifndef __LOGGER_HPP
#define __LOGGER_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <functional>

namespace logger
{

// log levels
enum level
{
    error = 0,
    warning,
    info,
    debug,
    debug1,
    debug2
};

class Logger
{
public:
    // Constructor: initialize buffer.
    Logger(std::string const & facility, level msg_level,
           std::string const & file_name, unsigned line_num, std::string const & func_name,
           std::ostream & os = std::clog)
        : _os_p(&os)
    {
        _oss << "= " << facility << "." << int(msg_level)
             << " " << file_name << ":" << line_num << " " << func_name << " ";
        _on_destruct = [&] () {
            _os_p->write(_oss.str().c_str(), _oss.str().size());
        };
    }
    // Constructor for exiting
    Logger(int exit_code,
           std::string const & file_name, unsigned line_num, std::string const & func_name,
           std::ostream & os = std::cerr)
        : _os_p(&os), _exit_code(exit_code)
    {
        _oss << file_name << ":" << line_num << " " << func_name << " ";
        _on_destruct = [&] () {
            _os_p->write(_oss.str().c_str(), _oss.str().size());
            if (_exit_code < 0)
            {
                std::abort();
            }
            else
            {
                std::exit(_exit_code);
            }
        };
    }
    // Constructor for throwing exceptions
    // first argument is only used to deduce the template argument type
    template <typename Exception>
    Logger(Exception const &,
           std::string const & file_name, unsigned line_num, std::string const & func_name,
           typename std::enable_if<std::is_base_of<std::exception, Exception>::value>::type * = 0)
    {
        _oss << file_name << ":" << line_num << " " << func_name << " ";
        _on_destruct = [&] () {
            throw Exception(_oss.str());
        };
    }
    // Destructor: dump buffer to output.
    ~Logger() noexcept(false)
    {
        _on_destruct();
    }
    // Produce l-value for output chaining.
    std::ostream & l_value() { return _oss; }

    // static methods for setting and getting facility log levels.
    static level get_default_level()
    {
        return default_level();
    }
    static void set_default_level(level l)
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lg(m);
        default_level() = l;
    }
    static void set_default_level(int l)
    {
        set_default_level(get_level(l));
    }
    static void set_default_level(std::string const & s)
    {
        set_default_level(get_level(s));
    }
    static level get_facility_level(std::string const & facility)
    {
        return (facility_level_map().count(facility) > 0?
                facility_level_map().at(facility) : get_default_level());
    }
    static void set_facility_level(std::string const & facility, level l)
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lg(m);
        facility_level_map()[facility] = l;
    }
    static void set_facility_level(std::string const & facility, int l)
    {
        set_facility_level(facility, get_level(l));
    }
    static void set_facility_level(std::string const & facility, std::string const & s)
    {
        set_facility_level(facility, get_level(s));
    }
    // static methods for setting log levels from command-line options
    static void set_level_from_option(std::string const & l, std::ostream * os_p = nullptr)
    {
        size_t i = l.find(':');
        if (i == std::string::npos)
        {
            set_default_level(l);
            if (os_p)
            {
                (*os_p) << "set default log level to: "
                        << static_cast<int>(Logger::get_default_level()) << std::endl;
            }
        }
        else
        {
            set_facility_level(l.substr(0, i), l.substr(i + 1));
            if (os_p)
            {
                (*os_p) << "set log level of '" << l.substr(0, i) << "' to: "
                        << static_cast<int>(Logger::get_facility_level(l.substr(0, i))) << std::endl;
            }
        }
    }
    static void set_levels_from_options(std::vector<std::string> const & v, std::ostream * os_p = nullptr)
    {
        for (auto const & l : v)
        {
            set_level_from_option(l, os_p);
        }
    }
    // public static utility functions (used by LOG macro)
    static level get_level(level l) { return l; }
    static level get_level(int i) { return static_cast<level>(i); }
    static level get_level(std::string const & s) { return level_from_string(s); }
    // public static member (used by LOG macro)
    static level& thread_local_last_level()
    {
        static thread_local level _last_level = error;
        return _last_level;
    }
private:
    std::ostringstream _oss;
    std::function<void()> _on_destruct;
    std::ostream * _os_p;
    int _exit_code;

    // private static data members
    static level & default_level()
    {
        static level _default_level = error;
        return _default_level;
    }
    static std::map<std::string, level> & facility_level_map()
    {
        static std::map<std::string, level> _facility_level_map;
        return _facility_level_map;
    }
    // private static utility functions
    static level level_from_string(std::string const & s)
    {
        std::istringstream iss(s + "\n");
        int tmp_int = -1;
        iss >> tmp_int;
        if (iss.good())
        {
            return level(tmp_int);
        }
        else
        {
            if (s == "error") return logger::error;
            else if (s == "warning") return logger::warning;
            else if (s == "info") return logger::info;
            else if (s == "debug") return logger::debug;
            else if (s == "debug1") return logger::debug1;
            else if (s == "debug2") return logger::debug2;
            else
            {
                std::ostringstream oss;
                oss << "could not parse log level: " << s;
                throw std::invalid_argument(oss.str());
            }
        }
    }
}; // class Logger

} //namespace logger

#define __FILENAME__ (std::string(__FILE__).find('/') != std::string::npos? std::string(__FILE__).substr(std::string(__FILE__).rfind('/') + 1) : std::string(__FILE__))

/**
 * LOG macro
 *
 * Synopsis:
 *   LOG(facility, level_spec, sink) << message
 *   LOG(facility, level_spec) << message
 *   LOG(level_spec) << message
 *
 *   `facility`   : string
 *   `level_spec` : integer, string, or logger level
 *   `sink`       : sink ostream
 * 
 * Log to `facility` at logger level `level_spec` and dump output to `sink`.
 * If sink is omitted, it defaults to std::clog.
 * If `facility` is omitted (logger has single argument), the macro LOG_FACILITY
 * is used instead, defaulting to "main".
 */

#define __LOG_3(facility, level_spec, sink)                                   \
    { using namespace logger; logger::Logger::thread_local_last_level() = logger::Logger::get_level(level_spec); } \
    if (logger::Logger::thread_local_last_level() > logger::Logger::get_facility_level(facility)) ; \
    else logger::Logger(facility, logger::Logger::thread_local_last_level(), __FILENAME__, __LINE__, __func__, sink).l_value()

#define __LOG_2(facility, level_spec)                                   \
    { using namespace logger; logger::Logger::thread_local_last_level() = logger::Logger::get_level(level_spec); } \
    if (logger::Logger::thread_local_last_level() > logger::Logger::get_facility_level(facility)) ; \
    else logger::Logger(facility, logger::Logger::thread_local_last_level(), __FILENAME__, __LINE__, __func__).l_value()

#define __LOG_1(level_spec) \
    __LOG_2(LOG_FACILITY, level_spec)

// we need 2-level indirection in order to trigger expansion after token pasting
// http://stackoverflow.com/questions/1597007/creating-c-macro-with-and-line-token-concatenation-with-positioning-macr
// http://stackoverflow.com/a/11763196/717706
#ifdef WIN32
#define __EXPAND(...) __VA_ARGS__
#define __LOG_aux2(N, ...) __EXPAND(__LOG_ ## N (__VA_ARGS__))
#else
#define __LOG_aux2(N, ...) __LOG_ ## N (__VA_ARGS__)
#endif

#define __LOG_aux1(N, ...) __LOG_aux2(N, __VA_ARGS__)

#define __NARGS_AUX(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, ...) _9

#ifdef WIN32
#define __NARGS(...) __EXPAND(__NARGS_AUX(__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0))
#else
#define __NARGS(...) __NARGS_AUX(__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0)
#endif

#ifndef LOG_FACILITY
#define LOG_FACILITY "main"
#endif

#define LOG(...) __LOG_aux1(__NARGS(__VA_ARGS__), __VA_ARGS__)

#define LOG_EXIT_(exit_code) logger::Logger((exit_code), __FILENAME__, __LINE__, __func__).l_value()
#define LOG_ABORT LOG_EXIT_(-1)
#define LOG_EXIT LOG_EXIT_(EXIT_FAILURE)

#define LOG_THROW_(Exception) logger::Logger(Exception(""), __FILENAME__, __LINE__, __func__).l_value()
#define LOG_THROW LOG_THROW_(std::runtime_error)

#endif

#ifdef SAMPLE_LOGGER

/*

Compile:

g++ -std=c++11 -D SAMPLE_LOGGER -x c++ logger.hpp -o sample-logger

Run:
./sample-logger info
./sample-logger info alt:debug1

*/

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "Use: " << argv[0] << " <log_level_setting> ..." << endl
             << "The program sends 5 log messages with decreasing priority (0=highest, 4=lowest)" << endl
             << "to 2 facilities \"main\" and \"alt\". Command-line arguments are interpreted as" << endl
             << "log facility level settings in the form [<facility>:]<level>." << endl;
        return EXIT_FAILURE;
    }
    for (int i = 1; i < argc; ++i)
    {
        cerr << "processing argument [" << argv[i] << "]" << endl;
        logger::Logger::set_level_from_option(argv[i], &cerr);
    }
    vector<string> const level_name{ "error", "warning", "info", "debug", "debug1", "debug2" };
    for (int l = 0; l < 5; ++l)
    {
        LOG(level_name[l]) << "message at level " << l << " (" << level_name[l]
                           << ") for facility main" << endl;
        LOG("alt", l) << "message at level " << l << " (" << level_name[l]
                      << ") for facility alt" << endl;
    }
}

#endif
