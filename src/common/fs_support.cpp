//---------------------------------------------------------
// Copyright 2016 Ontario Institute for Cancer Research
// Written by Matei David (matei.david@oicr.on.ca)
//---------------------------------------------------------
//
#include "fs_support.hpp"
#include <sys/stat.h>
#include <dirent.h>

//
bool is_directory(const std::string& file_name)
{
    auto dir = opendir(file_name.c_str());
    if(not dir) {
        return false;
    }
    closedir(dir);
    return true;
}

//
std::vector< std::string > list_directory(const std::string& file_name)
{
    std::vector< std::string > res;
    DIR* dir;
    struct dirent *ent;

    dir = opendir(file_name.c_str());
    if(not dir) {
        return res;
    }
    while((ent = readdir(dir)) != nullptr) {
        res.push_back(ent->d_name);
    }
    closedir(dir);
    return res;
}
