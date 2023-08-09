#pragma once

#include <iostream>
#include <string>
#include <vector>

// Includes ISL affine list/piecewise functions.
#include <isl/aff.h>
#include <isl/ilp.h>
// Includes ISL maps/binary relations.
#include <isl/map.h>
// Includes ISL ids and dspaces.
#include <isl/id.h>
#include <isl/space.h>
// Imports ISL val.
#include <isl/val.h>

long analyze_jumps(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_map *dist_func);
long analyze_jumps(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);
long analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_map *dist_func);
long analyze_latency(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);

// Defines debug variables from environment variables.
#include <string.h>
bool islIntermediates = (getenv("ISL_INTERMEDIATES") != NULL) &&
                        (strcmp(getenv("ISL_INTERMEDIATES"), "0") != 0);

// Defines a function to programatically generate an n-dimensional Manhattan distance function.
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims);
// Defines a function to programatically generate an n-circumference ring distance function.
std::string n_long_ring_metric(long ring_circumference);

// Defines debug dump function.
void dump(std::string str, isl_map *map)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_map_dump(map);
    }
}

void dump(std::string str, isl_pw_aff *pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_pw_aff_dump(pw_aff);
    }
}

void dump(std::string str, isl_multi_pw_aff *multi_pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_pw_aff_dump(multi_pw_aff);
    }
}

void dump(std::string str, isl_multi_val *multi_val)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_val_dump(multi_val);
    }
}