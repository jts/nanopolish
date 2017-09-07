//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Matei David (matei.david@oicr.on.ca)
//---------------------------------------------------------
//
#ifndef __FS_SUPPORT_HPP
#define __FS_SUPPORT_HPP

#include <string>
#include <vector>

// ref: http://stackoverflow.com/a/612176/717706
// return true if the given name is a directory
bool is_directory(const std::string& file_name);

// return a vector of files within the given directory
std::vector< std::string > list_directory(const std::string& file_name);

#endif
