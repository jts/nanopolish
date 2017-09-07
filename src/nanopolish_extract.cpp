//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Matei David (matei.david@oicr.on.ca)
//---------------------------------------------------------
//

//
// Getopt
//
#define SUBPROGRAM "extract"
#define LOG_FACILITY SUBPROGRAM

#include <iostream>
#include <sstream>
#include <getopt.h>

#include "nanopolish_read_db.h"
#include "nanopolish_extract.h"
#include "nanopolish_common.h"
#include "fs_support.hpp"
#include "logger.hpp"
#include "alg.hpp"

static const char *EXTRACT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2015 Ontario Institute for Cancer Research\n";

static const char *EXTRACT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTIONS] <fast5|dir>...\n"
"Extract reads in fasta format\n"
"\n"
"      --help                           display this help and exit\n"
"      --version                        display version\n"
"  -v, --verbose                        display verbose output\n"
"  -r, --recurse                        recurse into subdirectories\n"
"  -q, --fastq                          extract fastq (default: fasta)\n"
"  -t, --type=TYPE                      read type: template, complement, 2d, 2d-or-template, any\n"
"                                         (default: 2d-or-template)\n"
"  -b, --basecaller=NAME[:VERSION]      consider only data produced by basecaller NAME,\n"
"                                         optionally with given exact VERSION\n"
"  -o, --output=FILE                    write output to FILE\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose = 0;
    static bool recurse = false;
    static std::string read_type = "2d-or-template";
    static std::string basecaller_name;
    static std::string basecaller_version;
    static bool fastq = false;
    static std::string output_file;
    static std::deque< std::string > paths;
    static unsigned total_files_count = 0;
    static unsigned total_files_used_count = 0;
}
static std::ostream* os_p;

std::vector< std::pair< unsigned, std::string > >
get_preferred_basecall_groups(const fast5::File& f)
{
    bool have_2d = false;
    std::vector< std::pair< unsigned, std::string > > res;
    // check 2d
    if (opt::read_type == "any"
        or opt::read_type == "2d"
        or opt::read_type == "2d-or-template")
    {
        const auto& gr_l = f.get_basecall_strand_group_list(2);
        for (const auto& gr : gr_l)
        {
            auto bcd = f.get_basecall_group_description(gr);
            // Never use minknow basecalls
            if(bcd.name == "minknow") {
                continue;
            }

            if (not opt::basecaller_name.empty())
            {
                if ((bcd.name != opt::basecaller_name) or
                    (not opt::basecaller_version.empty() and
                     bcd.version != opt::basecaller_version))
                {
                    continue;
                }
            }

            if (f.have_basecall_fastq(2, gr)
                and f.have_basecall_events(0, gr)
                and f.have_basecall_events(1, gr))
            {
                have_2d = true;
                res.push_back(std::make_pair(2, gr));
                if (opt::read_type != "any")
                {
                    break;
                }
            }
        }
    }
    // check 1d
    for (unsigned st = 0; st < 2; ++st)
    {
        if (opt::read_type == "any"
            or (st == 0
                and (opt::read_type == "template"
                     or (not have_2d and opt::read_type == "2d-or-template")))
            or (st == 1
                and opt::read_type == "complement"))
        {
            const auto& gr_l = f.get_basecall_strand_group_list(st);
            for (const auto& gr : gr_l)
            {
                auto bcd = f.get_basecall_group_description(gr);
                // Never use minknow basecalls
                if(bcd.name == "minknow") {
                    continue;
                }

                if (not opt::basecaller_name.empty())
                {
                    if ((bcd.name != opt::basecaller_name) or
                        (not opt::basecaller_version.empty() and
                         bcd.version != opt::basecaller_version))
                    {
                        continue;
                    }
                }
                if (f.have_basecall_fastq(st, gr)
                    and f.have_basecall_events(st, gr))
                {
                    res.push_back(std::make_pair(st, gr));
                    if (opt::read_type != "any")
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
    LOG(debug) << fn << "\n";
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
            ++opt::total_files_count;
            // get preferred basecall groups
            auto l = get_preferred_basecall_groups(f);
            if (l.empty())
            {
                LOG(info) << "file [" << fn << "]: no basecalling data suitable for nanoplish\n";
                return;
            }
            ++opt::total_files_used_count;
            const auto& p = l.front();
            // get and parse fastq
            auto fq = f.get_basecall_fastq(p.first, p.second);
            auto fq_a = f.split_fq(fq);
            // construct name
            std::string name;
            std::istringstream(fq_a[0]) >> name;
            std::replace(name.begin(), name.end(), ':', '_');
            name += ":" + p.second + ":";
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
            if (not opt::fastq)
            {
                (*os_p)
                    << ">" << name << " " << base_fn << " " << fn << "\n"
                    << fq_a[1] << "\n";
            }
            else
            {
                (*os_p)
                    << "@" << name << " " << base_fn << " " << fn << "\n"
                    << fq_a[1] << "\n"
                    << "+" << fq_a[2] << "\n"
                    << fq_a[3] << "\n";
            }
        }
        catch (hdf5_tools::Exception& e)
        {
            LOG(warning) << fn << ": HDF5 error: " << e.what() << "\n";
        }
    } while (false);
} // process_file

void process_path(const std::string& path)
{
    LOG(info) << path << "\n";
    if (is_directory(path))
    {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list)
        {
            if (fn == "." or fn == "..") continue;
            std::string full_fn = path + "/" + fn;
            if (is_directory(full_fn))
            {
                if (opt::recurse)
                {
                    opt::paths.push_back(full_fn);
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
    }
    else // not is_directory; must be a fast5 file
    {
        if (fast5::File::is_valid_file(path))
        {
            process_file(path);
        }
        else
        {
            LOG(error) << "path [" << path << "] is neither a fast5 file nor a directory\n";
            exit(EXIT_FAILURE);
        }
    }
} // process_path

static const char* shortopts = "vrqt:o:b:";

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
    { "fastq",              no_argument,       NULL, 'q' },
    { "type",               required_argument, NULL, 't' },
    { "output",             required_argument, NULL, 'o' },
    { "basecaller",         required_argument, NULL, 'b' },
    { NULL, 0, NULL, 0 }
};

void parse_extract_options(int argc, char** argv)
{
    bool die = false;
    std::vector< std::string> log_level;
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
                log_level.push_back(arg.str());
                break;
            case 'v': opt::verbose++; break;
            case 'r': opt::recurse = true; break;
            case 'q': opt::fastq = true; break;
            case 't': arg >> opt::read_type; break;
            case 'o': arg >> opt::output_file; break;
        case 'b':
            arg >> opt::basecaller_name;
            auto i = opt::basecaller_name.find(':');
            if (i != std::string::npos)
            {
                opt::basecaller_version = opt::basecaller_name.substr(i + 1);
                opt::basecaller_name.resize(i);
            }
            break;
        }
    }
    // set log levels
    auto default_level = (int)logger::warning + opt::verbose;
    logger::Logger::set_default_level(default_level);
    logger::Logger::set_levels_from_options(log_level, &std::clog);

    if(opt::output_file.empty()) {
        std::cerr << SUBPROGRAM ": an output file (-o) is required\n";
        die = true;
    }

    // parse paths to process
    while (optind < argc)
    {
        std::string path = argv[optind++];
        while (path.size() > 1 and path[path.size() - 1] == '/')
        {
            path.resize(path.size() - 1);
        }
        opt::paths.push_back(path);
    }

    if (opt::paths.empty())
    {
        std::cerr << SUBPROGRAM ": no paths to process\n";
        die = true;
    }

    // check read type
    if (not (opt::read_type == "template"
             or opt::read_type == "complement"
             or opt::read_type == "2d"
             or opt::read_type == "2d-or-template"
             or opt::read_type == "any"))
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

    LOG(info) << "paths: " << alg::os_join(opt::paths, " ") << "\n";
    LOG(info) << "recurse: " << (opt::recurse? "yes" : "no") << "\n";
    LOG(info) << "read_type: " << opt::read_type << "\n";
    LOG(info) << "basecaller_name: " << opt::basecaller_name << "\n";
    LOG(info) << "basecaller_version: " << opt::basecaller_version << "\n";
}

int extract_main(int argc, char** argv)
{
    parse_extract_options(argc, argv);

    // Iterate over fast5 collection extracting the sequence reads
    // We do this in a block so the file is automatically closed
    // when ofs goes out of scope.
    {
        std::ofstream ofs;
        ofs.open(opt::output_file);
        os_p = &ofs;

        for (unsigned i = 0; i < opt::paths.size(); ++i)
        {
            process_path(opt::paths[i]);
        }
    }

    // Build the ReadDB from the output file
    ReadDB read_db;
    read_db.build(opt::output_file);
    read_db.save();

    std::clog << "[extract] found " << opt::total_files_count
              << " files, extracted " << opt::total_files_used_count
              << " reads\n";
    return 0;
}
