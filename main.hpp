#include <iostream>
#include <string>

#include <isl/map.h>

std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill);