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
    /* Defines an ISL distance function that calculates the manhattan distance
     * between two points. */
    isl_pw_aff *manhattan_metric = isl_pw_aff_read_from_str(
        p_ctx,
        "{"
            "[[xd, yd] -> [xs, ys]] -> [(xd - xs) + (yd - ys)] : "
                "xd >= xs and yd >= ys;"
            "[[xd, yd] -> [xs, ys]] -> [-(xd - xs) + -(yd - ys)] : "
                "xd < xs and yd < ys;"
            "[[xd, yd] -> [xs, ys]] -> [-(xd - xs) + (yd - ys)] : "
                "xd < xs and yd >= ys;"
            "[[xd, yd] -> [xs, ys]] -> [(xd - xs) + -(yd - ys)] : "
                "xd >= xs and yd < ys"
        "}"
    );
    std::cout << "manhattan_metric: " << std::endl;
    isl_pw_aff_dump(manhattan_metric);
    
    // Turns the above piecewise affine into a map.
    isl_map *manhattan_metric_map = isl_map_from_pw_aff(
        isl_pw_aff_copy(manhattan_metric)
    );
    std::cout << "manhattan_metric_map: " << std::endl;
    isl_map_dump(manhattan_metric_map);

  
    std::string result = analyze_latency(p_src_occupancy, p_dst_fill, manhattan_metric_map);
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
 * @param dist_func         The distance function to use, as a map.
 */
std::string analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_map *dist_func)
{
    // Prints out the inputs.
    std::cout << "p_src_occupancy: " << std::endl;
    isl_map_dump(p_src_occupancy);
    std::cout << "p_dst_fill: " << std::endl;
    isl_map_dump(p_dst_fill);
    std::cout << "dist_func: " << std::endl;
    isl_map_dump(dist_func);

    /* Invertsdst_fill such that data implies dst.
     * i.e. {[xd, yd] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *dst_fill_inverted = isl_map_reverse(
        isl_map_copy(p_dst_fill)
    );
    std::cout << "dst_fill_inverted: " << std::endl;
    isl_map_dump(dst_fill_inverted);
    /* Inverts src_occupancy such that data implies source.
     * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(
        isl_map_copy(p_src_occupancy)
    );
    std::cout << "src_occupancy_inverted: " << std::endl;
    isl_map_dump(src_occupancy_inverted);

    /* Takes the factored range of src_occupancy_inverted and dst_fill_inverted
     * to get {[d0, d1] -> [[xd, yd] -> [xs, ys]]} */
    isl_map *data_TO_dst_to_src = isl_map_range_product(
            isl_map_copy(dst_fill_inverted),
            isl_map_copy(src_occupancy_inverted)
    );
    std::cout << "data_TO_dst_to_src: " << std::endl;
    isl_map_dump(data_TO_dst_to_src);
    /* Wraps dst fill such that the binary relation implies data.
     * i.e. {[[xd, yd] -> [d0, d1]] -> [d0, d1]} */
    isl_map *dst_fill_wrapped = isl_map_range_map(
        isl_map_copy(p_dst_fill)
    );
    std::cout << "dst_fill_wrapped: " << std::endl;
    isl_map_dump(dst_fill_wrapped);

    /* Composites dst_fill_wrapped and data_to_dst_to_src to get
     * {[[xd, yd] -> [d0, d1]] -> [[xd, yd] -> [xs, ys]]} */
    isl_map *dst_to_data_TO_dst_to_src = isl_map_apply_range(
        isl_map_copy(dst_fill_wrapped),
        isl_map_copy(data_TO_dst_to_src)
    );
    std::cout << "dst_to_data_TO_dst_to_src: " << std::endl;
    isl_map_dump(dst_to_data_TO_dst_to_src);

    // Calculates the Manhattan distance for each [xd, yd] -> [d0, d1]
    isl_map *manhattan_distance = isl_map_apply_range(
        isl_map_copy(dst_to_data_TO_dst_to_src),
        isl_map_copy(dist_func)
    );
    
    std::cout << "manhattan_distance: " << std::endl;
    isl_map_dump(manhattan_distance);

    // Gets the minimum distance for each [xd, yd] -> [d0, d1]

    return "Coding In Progress...";
}