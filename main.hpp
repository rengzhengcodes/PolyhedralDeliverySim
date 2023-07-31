#include <iostream>
#include <string>

// Includes ISL affine list/piecewise functions.
#include <isl/aff.h>
// Includes ISL maps/binary relations.
#include <isl/map.h>

std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill);