#include <iostream>
#include <string>

// Includes ISL affine list/piecewise functions.
#include <isl/aff.h>
#include <isl/ilp.h>
// Includes ISL maps/binary relations.
#include <isl/map.h>

std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_map *dist_func);

// Defines debug variables from environment variables.
#include <string.h>
bool islIntermediates = (getenv("ISL_INTERMEDIATES") != NULL) &&
                        (strcmp(getenv("ISL_INTERMEDIATES"), "0") != 0);

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