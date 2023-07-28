#include "main.hpp"

int main(int argc, char* argv[])
{
    // Allocates space to the isl context and returns a pointer to it.
    isl_ctx *p_ctx = isl_ctx_alloc();
    /* Reads a map from a string relating source location and the data occupied
     * within. The map is represented using a binary relational diagram where
     * source implies data.*/
    isl_map *p_src_occupancy = isl_map_read_from_str(
        p_ctx,
        "{ [xs, ys] -> [d0, d1] : d0=xs and 0 <= d1 < 8 and 0 <= xs < 8 and 0 <= ys < 8 }"
    );
    /* Reads a map from a string relating destination location and the data
     * requested. The map is represented using a binary relational diagram where
     * source implies data. */
    isl_map *p_dst_fill = isl_map_read_from_str(
        p_ctx,
        "{ [xd, yd] -> [d0, d1] : d0=xd and d1=yd and 0 <= xd < 8 and 0 <= yd < 8 }"
    );
  
    std::string result = analyze_latency(p_src_occupancy, p_dst_fill);
    std::cout << result << std::endl;
}

/**
 * 
 */
std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill)
{
    return std::string("NOT IMPLEMENTED YET!!!");
}