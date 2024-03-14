#include "latency.hpp"

/// NOTES FOR NON-TREE MULTICAST SCENARIO
// - Load balancing issues for multiple minimally distant sources.
// - Compose minimal distances with the other set to remove non-minimal pairs then
// move on with the rest of the algorithm.
__isl_give isl_pw_qpolynomial *identify_minimal_pairs( 
    __isl_take isl_map *src_occupancy, 
    __isl_take isl_map *dst_fill, 
    __isl_take isl_map *dist_func
) {
    /* Makes [[dst -> data] -> dst] -> [data] */
    isl_set* wrapped_dst_fill = isl_map_wrap(dst_fill);
    isl_map* wrapped_fill_identity =
        isl_map_identity(isl_space_map_from_set(isl_set_get_space(
            wrapped_dst_fill
        )));
    dump("wrapped_fill_identity", wrapped_fill_identity);
    wrapped_fill_identity = isl_map_intersect_domain(
        wrapped_fill_identity,
        wrapped_dst_fill
    );
    dump("wrapped_fill_identity", wrapped_fill_identity);
    /* Makes [dst -> data] -> [dst -> data] */
    isl_map* uncurried_fill_identity = isl_map_uncurry(wrapped_fill_identity);
    dump("uncurried_fill_identity", uncurried_fill_identity);
    /* Inverts src_occupancy such that data implies source.
    * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(src_occupancy);

    isl_map* dst_to_data_to_dst_TO_src = isl_map_apply_range(
        uncurried_fill_identity,
        src_occupancy_inverted
    );
    dump("dst_to_data_to_dst_TO_src", dst_to_data_to_dst_TO_src);

    isl_map* dst_to_data_TO_dst_to_src =
        isl_map_curry(dst_to_data_to_dst_TO_src);

    // Calculates the distance of all the dst-src pairs with matching data.
    isl_map *distances_map = isl_map_apply_range(
        dst_to_data_TO_dst_to_src, dist_func
    );
    dump("distances_map", distances_map);

    // Converts the distances map to a piecewise affine.
    isl_map *lexmin_distances = isl_map_lexmin(distances_map);
    isl_multi_pw_aff *dirty_distances_aff =isl_multi_pw_aff_from_pw_multi_aff(isl_pw_multi_aff_from_map(lexmin_distances));
    // std::cout << "dirty_distances_aff" <<  isl_multi_pw_aff_to_str(dirty_distances_aff) << std::endl;
    assert(isl_multi_pw_aff_size(dirty_distances_aff) == 1);
    isl_pw_aff *distances_aff = isl_multi_pw_aff_get_at(dirty_distances_aff, 0);
    isl_multi_pw_aff_free(dirty_distances_aff);

    // Converts to a pw_qpolynomial for easier processing later.
    isl_pw_qpolynomial *distances_pwqp = isl_pw_qpolynomial_from_pw_aff(distances_aff);
    return distances_pwqp;
}   

isl_pw_qpolynomial *identify_minimal_pairs(
    const std::string& src_occupancy, 
    const std::string& dst_fill, 
    const std::string& dist_func
) {
    // Creates a new isl context.
    isl_ctx *p_ctx = isl_ctx_alloc();

    // Reads the string representations of the maps into isl objects.
    isl_map *p_src_occupancy = isl_map_read_from_str(
        p_ctx,
        src_occupancy.c_str()
    );
    isl_map *p_dst_fill = isl_map_read_from_str(
        p_ctx,
        dst_fill.c_str()
    );
    isl_map *p_dist_func = isl_map_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // Calls the isl version of analyze_latency.
    isl_pw_qpolynomial *ret = identify_minimal_pairs(
        p_src_occupancy,
        p_dst_fill,
        p_dist_func
    );

    // Frees the isl objects.
    isl_ctx_free(p_ctx);

    return ret;
}

int main(int argc, char* argv[])
{
    int M_int = 1024;
    int N_int = 1024;
    std::string M = std::to_string(M_int);
    std::string N = std::to_string(N_int);
    std::vector<int> D_vals({1});//, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024});
    clock_t start, end;
    double cpu_time_used;
    for (int D_int : D_vals) {
        start = clock();
        std::string D = std::to_string(D_int);
        // Defines the src occupancy map as a string.
        std::string src_occupancy = "{[xs, ys] -> [a, b] : ("+D+"*xs)%"+M+" <= a <= ("+
                                    D+"*xs+"+D+"-1)%"+M+" and b=ys and 0 <= xs < "+M+
                                    " and 0 <= ys < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";
        // Defines the dst fill map as a string.
        std::string dst_fill =  "{[xd, yd] -> [a, b] : b=yd and 0 <= xd < "+M+
                                " and 0 <= yd < "+N+" and 0 <= a < "+M+" and 0 <= b < "+N+" }";

        // Defines the distance function string.
        std::string dist_func_str = nd_manhattan_metric({"xs", "ys"}, {"xd", "yd"});
    
        identify_minimal_pairs(src_occupancy, dst_fill, dist_func_str);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    }
}