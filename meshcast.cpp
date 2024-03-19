#include "latency.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

/// NOTES FOR NON-TREE MULTICAST SCENARIO
// - Load balancing issues for multiple minimally distant sources.
// - Compose minimal distances with the other set to remove non-minimal pairs then
// move on with the rest of the algorithm.
__isl_give isl_map *identify_mesh_casts( 
    __isl_take isl_map *src_occupancy, 
    __isl_take isl_map *dst_fill, 
    __isl_take isl_map *dist_func
) {
    /* Makes [[dst -> data] -> dst] -> [data] */
    isl_set *wrapped_dst_fill = isl_map_wrap(dst_fill);
    isl_map *wrapped_fill_identity =
        isl_map_identity(isl_space_map_from_set(isl_set_get_space(
            wrapped_dst_fill
        )));
    DUMP(wrapped_fill_identity);
    wrapped_fill_identity = isl_map_intersect_domain(
        wrapped_fill_identity,
        wrapped_dst_fill
    );
    DUMP(wrapped_fill_identity);
    /* Makes [dst -> data] -> [dst -> data] */
    isl_map *uncurried_fill_identity = isl_map_uncurry(wrapped_fill_identity);
    DUMP(uncurried_fill_identity);
    /* Inverts src_occupancy such that data implies source.
    * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(src_occupancy);

    isl_map *dst_to_data_to_dst_TO_src = isl_map_apply_range(
        uncurried_fill_identity,
        src_occupancy_inverted
    );
    DUMP(dst_to_data_to_dst_TO_src);

    isl_map *dst_to_data_TO_dst_to_src =
        isl_map_curry(dst_to_data_to_dst_TO_src);

    // Calculates the distance of all the dst-src pairs with matching data.
    isl_map *distances_map = isl_map_apply_range(
        isl_map_range_map(dst_to_data_TO_dst_to_src), dist_func
    );
    DUMP(distances_map);

    // Gets the minimal distance pairs.
    isl_map *lexmin_distances = isl_map_lexmin(distances_map);
    // Isolates the relevant domain portions.
    isl_map *minimal_pairs = isl_set_unwrap(isl_map_domain(lexmin_distances));
    DUMP(minimal_pairs);
    // Isolates the multicast networks.
    isl_map *multicast_networks = isl_map_curry(minimal_pairs);
    multicast_networks = isl_set_unwrap(isl_map_range(multicast_networks));
    DUMP(multicast_networks);

    return minimal_pairs;
}   

isl_map *identify_mesh_casts(
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
    isl_map *ret = identify_mesh_casts(
        p_src_occupancy,
        p_dst_fill,
        p_dist_func
    );

    // Frees the isl objects.
    isl_ctx_free(p_ctx);

    return ret;
}

__isl_give isl_map *cost_mesh_cast(
    __isl_take isl_map *mesh_cast_networks,
    __isl_take isl_map *dist_func
) {
    DUMP(mesh_cast_networks);
    DUMP(dist_func);
    
    // Calculates the cost per pair per datum.
    isl_map *cost_per_pair_per_datum = isl_map_range_map(mesh_cast_networks);
    cost_per_pair_per_datum = isl_map_apply_range(cost_per_pair_per_datum, dist_func);
    DUMP(cost_per_pair_per_datum);

    return cost_per_pair_per_datum;
}

/**
 * @return The cost per datum of each network.
 */
isl_map *cost_mesh_cast(
    const std::string& mesh_cast_networks,
    const std::string& dist_func
) {
    // Creates a new isl context.
    isl_ctx *p_ctx = isl_ctx_alloc();

    // Reads the string representations of the maps into isl objects.
    isl_map *p_mesh_cast_networks = isl_map_read_from_str(
        p_ctx,
        mesh_cast_networks.c_str()
    );
    isl_map *p_dist_func = isl_map_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // Calls the isl version of analyze_latency.
    isl_map *ret = cost_mesh_cast(
        p_mesh_cast_networks,
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
        std::string dist_func_str = R"DIST({
            [[xd, yd] -> [xs, ys]] -> [(xd - xs) + (yd - ys)] : 
                xd >= xs and yd >= ys;
            [[xd, yd] -> [xs, ys]] -> [-(xd - xs) + -(yd - ys)] : 
                xd < xs and yd < ys;
            [[xd, yd] -> [xs, ys]] -> [-(xd - xs) + (yd - ys)] : 
                xd < xs and yd >= ys;
            [[xd, yd] -> [xs, ys]] -> [(xd - xs) + -(yd - ys)] : 
                xd >= xs and yd < ys
            })DIST";
    
        auto mcs = identify_mesh_casts(src_occupancy, dst_fill, dist_func_str);
        DUMP(mcs);
        auto res = cost_mesh_cast(isl_map_to_str(mcs), dist_func_str);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    }
}