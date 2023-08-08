#include "latency.hpp"

int main(int argc, char* argv[])
{
    // Defines the src occupancy map as a string.
    std::string src_occupancy = "{ [xs] -> [d0] : 0 <= xs < 8 and d0 = xs }";
    // Defines the dst fill map as a string.
    std::string dst_fill = "{ [xd] -> [d0] : 0 <= d0 < 8 and xd = 0 }";

    // Defines the torus circumference.
    int torus_circumference = 8;
    // Defines the distance function string.
    std::string dist_func_str = n_long_ring_metric(torus_circumference);

    std::cout << dist_func_str << std::endl;
  
    std::string result = analyze_latency(src_occupancy, dst_fill, dist_func_str);
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
std::string analyze_latency (
    isl_map *p_src_occupancy, 
    isl_map *p_dst_fill, 
    isl_map *dist_func
) {
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

    // Restricts the range such that xd' = xd and yd' = yd.
    if (islIntermediates) std::cout << '\n' << std::endl;
    for (int i = 0; i < isl_map_dim(p_dst_fill, isl_dim_in); i++)
    {
        /* Restricts the ith element of the output by equating it to the ith
         * element of the input. Treats input and output as if it were flat. */
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

    // Gets the string representation of the maximum minimum distance.
    std::string result = isl_multi_val_to_str(max_min_distance);

    // Frees the isl intermediate objects.
    isl_map_free(dst_fill_inverted);
    isl_map_free(src_occupancy_inverted);
    isl_map_free(data_TO_dst_to_src);
    isl_map_free(dst_fill_wrapped);
    isl_map_free(dst_to_data_TO_dst_to_src);
    isl_map_free(manhattan_distance);
    isl_multi_pw_aff_free(min_distance);
    isl_multi_val_free(max_min_distance);

    return result;
}

/**
 * A wrapper for analyze_latency that takes in strings instead of isl objects.
 * 
 * @param src_occupancy     A string representation of a map relating source
 *                         location and the data occupied.
 * @param dst_fill          A string representation of a map relating destination
 *                         location and the data requested.
 * @param dist_func         A string representation of a distance function to use.
 *
 * @return                  A string representation of the maximum latency.
 */
std::string analyze_latency (
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
    isl_pw_aff *p_dist_func = isl_pw_aff_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // Turns dist_func into a map.
    isl_map *p_dist_func_map = isl_map_from_pw_aff(
        isl_pw_aff_copy(p_dist_func)
    );

    // Dumps the isl objects.
    dump("\np_src_occupancy: ", p_src_occupancy);
    dump("p_dst_fill: ", p_dst_fill);
    dump("p_dist_func: ", p_dist_func);
    dump("p_dist_func_map: ", p_dist_func_map);

    // Frees p_dist_func as it's now no longer needed once the map is formed.
    isl_pw_aff_free(p_dist_func);

    // Calls the isl version of analyze_latency.
    std::string result = analyze_latency(
        p_src_occupancy,
        p_dst_fill,
        p_dist_func_map
    );

    // Frees the isl objects.
    isl_map_free(p_src_occupancy);
    isl_map_free(p_dst_fill);
    isl_map_free(p_dist_func_map);
    isl_ctx_free(p_ctx);

    return result;
}

/**
 * Defines the n-dimensional Manhattan distance function. This is done programatically
 * as ISL does not have an absolute value function.
 * 
 * @pre             src_dims.size() == dst_dims.size()
 * 
 * @param src_dims  A vector of strings representing the source dimensions.
 * @param dst_dims  A vector of strings representing the destination dimensions.
 * 
 * @return          A piecewise affine function string representing the Manhattan distance.
 */
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims)
{
    // Creates a new isl context.
    isl_ctx *p_ctx = isl_ctx_alloc();

    /* Creates isl_ids for the src and dst dimensions. This is to be used for 
     * the map as unique identifiers. */
    std::vector<isl_id*> src_ids;
    std::vector<isl_id*> dst_ids;
    for (int i = 0; i < src_dims.size(); i++)
    {
        src_ids.push_back(
            isl_id_alloc(
                p_ctx,
                src_dims[i].c_str(),
                nullptr
            )
        );
        dst_ids.push_back(
            isl_id_alloc(
                p_ctx,
                dst_dims[i].c_str(),
                nullptr
            )
        );
    }

    // Programmatically binds the dst and src dimensions to their respective spaces.
    isl_space *p_in_space = isl_space_alloc(
        p_ctx,
        dst_dims.size(),
        0,
        0
    );
    isl_space *p_out_space = isl_space_alloc(
        p_ctx,
        src_dims.size(),
        0,
        0
    );
    for (int i = 0; i < src_dims.size(); i++)
    {
        p_in_space = isl_space_set_dim_id(
            p_in_space,
            isl_dim_param,
            i,
            isl_id_copy(dst_ids[i])
        );
        p_out_space = isl_space_set_dim_id(
            p_out_space,
            isl_dim_param,
            i,
            isl_id_copy(src_ids[i])
        );
    }
    isl_space_dump(p_in_space);

    // Vector of all the affines for each part of the piecewise function.
    std::vector<isl_aff*> affs;
    // Constructs all the absolute value affines per dimension.
    for (int i = 0; i < src_dims.size(); i++)
    {
        // Constructs the affine for the src.
        isl_aff *s_aff = isl_aff_var_on_domain(
            isl_local_space_from_space(p_out_space),
            isl_dim_in,
            i
        );
        // Constructs the affine for the dst.
        isl_aff *d_aff = isl_aff_var_on_domain(
            isl_local_space_from_space(p_in_space),
            isl_dim_in,
            i
        );
        // Subtracts the dst. 
        isl_aff *p_aff = isl_aff_sub(s_aff, d_aff);
        isl_aff_dump(p_aff);
    }


    return "Manhattan Scaffolding...\n"
    "          __  __                                             \n"
    "         |. ||. |    .|                                      \n"
    "         || ||| |    | |                W                    \n"
    "         |: ||: |    |'|               [ ]         ._____    \n"
    "         |  ||  |   |  |     .--'|      3   .---\"| |.   |'   \n"
    "     _   |  ||  |-. |  | __  |.  |     /|  _|__  | ||   |__  \n"
    "  .-'|  _|  ||  | ||   '-  | ||    \\|// / |   |' | |    | |' \n"
    "  |' | |.|  ||  | ||       '-'    -( )-|  |   |  | |    | |  \n"
    "__|  '-' '  ''  ' ""       '       J V |  `   -  |_'    ' |__\n"
    "                             ___  '    /                     \n"
    "                             \\  \\/    |                      \n"
    "Hilsen, Peer W Hansen--                                      \n";
}

/**
 * Calculates the latency of a memory access on a ring.
 * 
 * @param n   The circumference of the torus. 
 */
std::string n_long_ring_metric(long n)
{
    // Creates a new isl context.
    isl_ctx *p_ctx = isl_ctx_alloc();

    /* Creates isl_ids for the src and dst dimensions. This is to be used for
     * the map as unique identifiers. */
    isl_id *src_id = isl_id_alloc(p_ctx, "src", NULL);
    isl_id *dst_id = isl_id_alloc(p_ctx, "dst", NULL);

    // Creates the space in which distance calculations are done.
    isl_space *p_dist_space = isl_space_set_alloc(p_ctx, 0, 2);
    
    // Sets the dst and src as members of the space.
    p_dist_space = isl_space_set_dim_id(isl_space_copy(p_dist_space), isl_dim_set, 0, dst_id);
    p_dist_space = isl_space_set_dim_id(isl_space_copy(p_dist_space), isl_dim_set, 1, src_id);

    // Creates p_dist as a local space.
    isl_local_space *p_dist_local = isl_local_space_from_space(p_dist_space);

    // Creates the src and dst affines.
    isl_pw_aff *src_aff = isl_pw_aff_var_on_domain(isl_local_space_copy(p_dist_local), isl_dim_set, 1);
    isl_pw_aff *dst_aff = isl_pw_aff_var_on_domain(isl_local_space_copy(p_dist_local), isl_dim_set, 0);

    // Subtracts the dst from the src.
    isl_pw_aff *src_sub_dst_aff = isl_pw_aff_sub(isl_pw_aff_copy(src_aff), isl_pw_aff_copy(dst_aff));
    // Subtracts the src from the dst.
    isl_pw_aff *dst_sub_src_aff = isl_pw_aff_sub(isl_pw_aff_copy(dst_aff), isl_pw_aff_copy(src_aff));

    // Creates an isl_val for the torus circumference.
    isl_val *p_circumference = isl_val_int_from_si(p_ctx, n);

    // Creates the moduli for both differences.
    isl_pw_aff *src_sub_dst_mod_n_aff = isl_pw_aff_mod_val(src_sub_dst_aff, isl_val_copy(p_circumference));
    isl_pw_aff *dst_sub_src_mod_n_aff = isl_pw_aff_mod_val(dst_sub_src_aff, isl_val_copy(p_circumference));

    // Combines the moduli into a single piecewise affine.
    isl_pw_aff *p_dist = isl_pw_aff_min(src_sub_dst_mod_n_aff, dst_sub_src_mod_n_aff);
    // Grabs the return value as a string.
    std::string ret = isl_pw_aff_to_str(p_dist);

    return ret;
}