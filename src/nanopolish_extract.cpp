//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Matei David (matei.david@oicr.on.ca)
//---------------------------------------------------------
//
#include "nanopolish_extract.h"

#include <iostream>
#include <sstream>
#include <getopt.h>

#include "nanopolish_common.h"
#include "fs_support.hpp"
#include "logger.hpp"
#include "alg.hpp"
#include <fast5.hpp>

//
// Getopt
//
#define SUBPROGRAM "extract"
#define LOG_FACILITY SUBPROGRAM

static const char *EXTRACT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *EXTRACT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] <dir>...\n"
"Extract reads in fasta format\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display version\n"
"  -v, --verbose                        display verbose output\n"
"  -r, --recurse                        recurse into subdirectories\n"
"  -t, --type=TYPE                      read type: template, complement, 2d, 2d-or-template, all\n"
"                                         (default: 2d-or-template)\n"
"  -o, --output=FILE                    write output to FILE (default: stdout)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose = 0;
    static bool recurse = false;
    static std::string read_type = "2d-or-template";
    static std::string output_file;
    static std::deque< std::string > dirs;
}
static std::ostream* os_p;

std::vector< std::pair< unsigned, std::string > >
get_preferred_basecall_groups(const fast5::File& f, const std::string& read_type)
{
    bool have_2d = false;
    std::vector< std::pair< unsigned, std::string > > res;
    // check 2d
    if (read_type == "all"
        or read_type == "2d"
        or read_type == "2d-or-template")
    {
        const auto& gr_l = f.get_basecall_strand_group_list(2);
        for (const auto& gr : gr_l)
        {
            if (f.have_basecall_fastq(2, gr)
                and f.have_basecall_events(0, gr)
                and f.have_basecall_events(1, gr))
            {
                have_2d = true;
                res.push_back(std::make_pair(2, gr));
                if (read_type != "all")
                {
                    break;
                }
            }
        }
    }
    // check 1d
    for (unsigned st = 0; st < 2; ++st)
    {
        if (opt::read_type == "all"
            or (st == 0
                and (opt::read_type == "template"
                     or (not have_2d and opt::read_type == "2d-or-template")))
            or (st == 1
                and opt::read_type == "complement"))
        {
            const auto& gr_l = f.get_basecall_strand_group_list(st);
            for (const auto& gr : gr_l)
            {
                if (f.have_basecall_fastq(st, gr)
                    and f.have_basecall_events(st, gr))
                {
                    res.push_back(std::make_pair(st, gr));
                    if (read_type != "all")
                    {
                        break;
                    }
                }
            }
        }
    }
    LOG(debug)
        << "preferred_groups: "
        << alg::os_join(res, ", ", [] (const decltype(res)::value_type& p) {
                std::ostringstream oss;
                oss << "(" << p.first << ":" << p.second << ")";
                return oss.str();
            })
        << "\n";
    return res;
} // get_preferred_basecall_groups

void process_file(const std::string& fn)
{
    LOG(debug) << "processing_file: " << fn << "\n";
    auto pos = fn.find_last_of('/');
    std::string base_fn = (pos != std::string::npos? fn.substr(pos + 1) : fn);
    if (base_fn.substr(base_fn.size() - 6) == ".fast5")
    {
        base_fn.resize(base_fn.size() - 6);
    }
    fast5::File f;
    do
    {
        try
        {
            // open file
            f.open(fn);
            // get preferred basecall groups
            auto l = get_preferred_basecall_groups(f, opt::read_type);
            if (l.empty()) return;
            const auto& p = l.front();
            // get and parse fastq
            auto fq = f.get_basecall_fastq(p.first, p.second);
            auto fq_a = f.split_fq(fq);
            // construct name
            auto pos = fq_a[0].find_first_of('_');
            std::string name = fq_a[0].substr(0, pos) + ":Basecall_" + p.second + ":";
            if (p.first == 0)
            {
                name += "template";
            }
            else if (p.first == 1)
            {
                name += "complement";
            }
            else
            {
                name += "2d";
            }
            (*os_p) << ">" << name << " " << base_fn << " " << fn << "\n"
                    << fq_a[1] << "\n";
        }
        catch (hdf5_tools::Exception& e)
        {
            LOG(warning) << fn << ": HDF5 error: " << e.what() << "\n";
        }
    } while (false);
} // process_file

void process_dir(const std::string& dir)
{
    LOG(info) << "processing_dir: " << dir << "\n";
    auto dir_list = list_directory(dir);
    for (const auto& fn : dir_list)
    {
        if (fn == "." or fn == "..") continue;
        std::string full_fn = dir + (dir[dir.size() - 1] != '/'? "/" : "") + fn;
        if (is_directory(full_fn))
        {
            if (opt::recurse)
            {
                opt::dirs.push_back(full_fn);
            }
            else
            {
                LOG(info) << "ignoring_subdir: " << full_fn << "\n";
            }
        }
        else if (fast5::File::is_valid_file(full_fn))
        {
            process_file(full_fn);
        }
        else
        {
            LOG(info) << "ignoring_file: " << full_fn << "\n";
        }
    }
} // process_dir

static const char* shortopts = "vrt:o:";

enum {
    OPT_HELP = 1,
    OPT_VERSION,
    OPT_LOG_LEVEL,
};

static const struct option longopts[] = {
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "log-level",          required_argument, NULL, OPT_LOG_LEVEL },
    { "verbose",            no_argument,       NULL, 'v' },
    { "recurse",            no_argument,       NULL, 'r' },
    { "type",               required_argument, NULL, 't' },
    { "output",             required_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

void parse_extract_options(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case OPT_HELP:
                std::cout << EXTRACT_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << EXTRACT_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_LOG_LEVEL:
                Logger::set_level_from_option(arg.str());
                break;
            case 'v': opt::verbose++; break;
            case 'r': opt::recurse = true; break;
            case 't': arg >> opt::read_type; break;
            case 'o': arg >> opt::output_file; break;
        }
    }
    // parse dirs to process
    while (optind < argc)
    {
        auto dir = argv[optind++];
        if (not is_directory(dir))
        {
            std::cerr << SUBPROGRAM ": not a directory: " << dir << "\n";
            die = true;
        }
        opt::dirs.push_back(dir);
    }
    if (opt::dirs.empty())
    {
        std::cerr << SUBPROGRAM ": no directories to process\n";
        die = true;
    }
    // check read type
    if (not (opt::read_type == "template"
             or opt::read_type == "complement"
             or opt::read_type == "2d"
             or opt::read_type == "2d-or-template"
             or opt::read_type == "all"))
    {
        std::cerr << SUBPROGRAM ": invalid read type: " << opt::read_type << "\n";
        die = true;
    }
    // die if errors
    if (die)
    {
        std::cerr << "\n" << EXTRACT_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    LOG(info) << "dirs: " << alg::os_join(opt::dirs, " ") << "\n";
    LOG(info) << "recurse: " << (opt::recurse? "yes" : "no") << "\n";
    LOG(info) << "read_type: " << opt::read_type << "\n";
}

int extract_main(int argc, char** argv)
{
    parse_extract_options(argc, argv);
    std::ofstream ofs;
    if (not opt::output_file.empty())
    {
        ofs.open(opt::output_file);
        os_p = &ofs;
    }
    else
    {
        os_p = &std::cout;
    }
    for (unsigned i = 0; i < opt::dirs.size(); ++i)
    {
        process_dir(opt::dirs[i]);
    }
    return 0;
}
