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
        "{ [xs, ys] -> [d0, d1] : d0=xs and d1=ys and 0 <= xs < 8 and 0 <= ys < 8 }"
    );
    /* Reads a map from a string relating destination location and the data
     * requested. The map is represented using a binary relational diagram where
     * source implies data. */
    isl_map *p_dst_fill = isl_map_read_from_str(
        p_ctx,
        "{ [xd, yd] -> [d0, d1] : d0=xd and 0 <= d1 < 8 and 0 <= xd < 8 and 0 <= yd < 8 }"
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
    dump("manhattan_metric: ", manhattan_metric);
    
    // Turns the above piecewise affine into a map.
    isl_map *manhattan_metric_map = isl_map_from_pw_aff(
        isl_pw_aff_copy(manhattan_metric)
    );
    dump("manhattan_metric_map: ", manhattan_metric_map);

  
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
    dump("\np_src_occupancy: ", p_src_occupancy);
    dump("p_dst_fill: ", p_dst_fill);
    dump("dist_func: ", dist_func);

    /* Invertsdst_fill such that data implies dst.
     * i.e. {[xd, yd] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *dst_fill_inverted = isl_map_reverse(
        isl_map_copy(p_dst_fill)
    );
    dump("\ndst_fill_inverted: ", dst_fill_inverted);
    /* Inverts src_occupancy such that data implies source.
     * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(
        isl_map_copy(p_src_occupancy)
    );
    dump("src_occupancy_inverted: ", src_occupancy_inverted);

    /* Takes the factored range of src_occupancy_inverted and dst_fill_inverted
     * to get {[d0, d1] -> [[xd, yd] -> [xs, ys]]} */
    isl_map *data_TO_dst_to_src = isl_map_range_product(
            isl_map_copy(dst_fill_inverted),
            isl_map_copy(src_occupancy_inverted)
    );
    dump("\ndata_TO_dst_to_src: ", data_TO_dst_to_src);
    /* Wraps dst fill such that the binary relation implies data.
     * i.e. {[[xd, yd] -> [d0, d1]] -> [d0, d1]} */
    isl_map *dst_fill_wrapped = isl_map_range_map(
        isl_map_copy(p_dst_fill)
    );
    dump("dst_fill_wrapped: ", dst_fill_wrapped);

    /* Composites dst_fill_wrapped and data_to_dst_to_src to get
     * {[[xd, yd] -> [d0, d1]] -> [[xd', yd'] -> [xs, ys]]} */
    isl_map *dst_to_data_TO_dst_to_src = isl_map_apply_range(
        isl_map_copy(dst_fill_wrapped),
        isl_map_copy(data_TO_dst_to_src)
    );
    dump("\ndst_to_data_TO_dst_to_src: ", dst_to_data_TO_dst_to_src);

    if (islIntermediates) std::cout << '\n' << std::endl;
    // Restricts the range such that xd' = xd and yd' = yd.
    for (int i = 0; i < isl_map_dim(p_dst_fill, isl_dim_in); i++)
    {
        dst_to_data_TO_dst_to_src = isl_map_equate(
            dst_to_data_TO_dst_to_src,
            isl_dim_in,
            i,
            isl_dim_out,
            i
        );
        dump("dst_to_data_TO_dst_to_src_restricting_"+std::to_string(i)+": ", 
             dst_to_data_TO_dst_to_src);
    };
    dump("\ndst_to_data_TO_dst_to_src_restricted: ", dst_to_data_TO_dst_to_src);

    /* Computes the manhattan distance between the destination for a data and
     * a source for that data. */
    isl_map *manhattan_distance = isl_map_apply_range(
        isl_map_copy(dst_to_data_TO_dst_to_src),
        isl_map_copy(dist_func)
    );
    dump("\nmanhattan_distance: ", manhattan_distance);

    // Computes the minimum distance from every source to every destination.
    isl_multi_pw_aff *min_distance = isl_map_min_multi_pw_aff(
        isl_map_copy(manhattan_distance)
    );
    dump("\nmin_distance: ", min_distance);

    // Computes the maximum of minimum distances for every data.
    isl_multi_val *max_min_distance = isl_multi_pw_aff_max_multi_val(
        isl_multi_pw_aff_copy(min_distance)
    );
    dump("\nmax_min_distance: ", max_min_distance);

    return "Coding In Progress...";
}