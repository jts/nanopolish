#ifndef NANOPOLISH_EXTRACT_H
#define NANOPOLISH_EXTRACT_H

#include <string>
#include <vector>

#include <fast5.hpp>

std::vector< std::pair< unsigned, std::string > >
get_preferred_basecall_groups(const fast5::File& f, const std::string& read_type = "2d-or-template");

int extract_main(int argc, char** argv);

#endif
