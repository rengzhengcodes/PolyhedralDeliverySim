#include "main.hpp"

int main(int argc, char* argv[])
{
    // Allocates space to the isl context and returns a pointer to it.
    isl_ctx *p_ctx = isl_ctx_alloc();
    /* Reads a map from a string relating source location and the data occupied
     * within. The map is represented using a binary relational diagram where
     * source implies data. */
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
 * Analyzes the latency of a memory access by finding the minimum path from
 * every source to every destination for a particular data, then taking the max
 * of all the minimum paths for that data, then taking the max of all the latencies
 * for all the data.
 * 
 * @param p_src_occupancy   A map relating source location and the data occupied.
 * @param p_dst_fill        A map relating destination location and the data 
 *                          requested.
 */
std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill)
{
    // Wraps dst fill such that the binary relation implies data.
    isl_map *dst_fill_wrapped = isl_map_range_map(p_dst_fill);
    // Prints out dst_fill_wrapped in the terminal.
    // isl_map_dump(dst_fill_wrapped);
    
    // Inverts src_occupancy such that data implies source.
    isl_map *src_occupancy_inverted = isl_map_reverse(p_src_occupancy);
    // Prints out src_occupancy_inverted in the terminal.
    // isl_map_dump(src_occupancy_inverted);

    // Composites dst_fill_wrapped and src_occupancy_inverted such that data
    // implies source.
    isl_map *dst_fill_wrapped_composed = isl_map_apply_range(
        dst_fill_wrapped,
        src_occupancy_inverted
    );
    // Prints out dst_fill_wrapped_composed in the terminal.
    isl_map_dump(dst_fill_wrapped_composed);
    

    return "Coding In Progress...";
}