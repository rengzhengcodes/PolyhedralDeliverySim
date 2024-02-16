#include "latency.hpp"

struct qpolynomial_from_fold_info
{
  isl_pw_qpolynomial** pp_pwqp;
  isl_set* domain;
};

isl_stat fold_accumulator(isl_qpolynomial* qp, void* pwqp_out)
{
  auto p_info = static_cast<qpolynomial_from_fold_info*>(pwqp_out);
  auto p_pwqp = isl_pw_qpolynomial_from_qpolynomial(qp);
  p_pwqp = isl_pw_qpolynomial_intersect_domain(p_pwqp,
                                               isl_set_copy(p_info->domain));
  if (*p_info->pp_pwqp)
  {
    *p_info->pp_pwqp = isl_pw_qpolynomial_add(p_pwqp, *p_info->pp_pwqp);
  }
  else
  {
    *p_info->pp_pwqp = p_pwqp;
  }
  return isl_stat_ok;
}

isl_bool
pw_fold_accumulator(isl_set* set, isl_qpolynomial_fold* fold, void* pwqp_out)
{
  qpolynomial_from_fold_info info {
    .pp_pwqp = static_cast<isl_pw_qpolynomial**>(pwqp_out),
    .domain = set
  };
  isl_qpolynomial_fold_foreach_qpolynomial(
    fold,
    fold_accumulator,
    &info
  );
  return isl_bool_true;
}

__isl_give isl_pw_qpolynomial*
gather_pw_qpolynomial_from_fold(__isl_take isl_pw_qpolynomial_fold* pwqpf)
{
  isl_pw_qpolynomial* p_pwqp = nullptr;
  isl_pw_qpolynomial_fold_every_piece(
    pwqpf,
    pw_fold_accumulator,
    &p_pwqp
  );
  return p_pwqp;
}

int main(int argc, char* argv[])
{
    int M_int = 4;
    int N_int = 4;
    std::string M = std::to_string(M_int);
    std::string N = std::to_string(N_int);
    std::vector<int> D_vals({1, 2, 4});
    for (int D_int : D_vals) {
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
    
        // long latency = analyze_latency(src_occupancy, dst_fill, dist_func_str);
        // std::cout << "latency: " << latency << std::endl;
        long jumps = analyze_jumps(src_occupancy, dst_fill, dist_func_str);
        std::cout << "D: " << D << " | jumps:\t " << jumps << std::endl;
    }
}

/**
 * Minimizes the distance between every dst and src per data.
 * 
 * @param __isl_take p_src_occupancy    A map relating source location and the
 *                                      data occupied.
 * @param __isl_take p_dst_fill         A map relating destination location and
 *                                      the data requested.
 * @param __isl_take dist_func          The distance function to use, as a map.
 */ 
__isl_give isl_pw_qpolynomial *minimize_jumps(
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
    isl_pw_multi_aff *dirty_distances_aff = isl_pw_multi_aff_from_map(distances_map);
    assert(isl_pw_multi_aff_n_piece(dirty_distances_aff) == 1);
    isl_pw_aff *distances_aff = isl_pw_multi_aff_get_at(dirty_distances_aff, 1);
    isl_pw_multi_aff_free(dirty_distances_aff);

    // Converts to a pw_qpolynomial for easier processing later.
    isl_pw_qpolynomial *distances_pwqp = isl_pw_qpolynomial_from_pw_aff(distances_aff);

    return distances_pwqp;
}

/**
 * Analyzes the total jumps that a process takes given the source, destination,
 * and distance function.
 * 
 * @param p_src_occupancy    A map relating source location and the
 *                           data occupied.
 * @param p_dst_fill         A map relating destination location and
 *                           the data requested.
 * @param dist_func          The distance function to use, as a map.
 */
long analyze_jumps(
    __isl_take isl_map *src_occ, 
    __isl_take isl_map *dst_fill,
    __isl_take isl_map *dist_func
) {
    // Prints out inputs for debugging.
    dump("src_occupancy: ", src_occ);
    dump("dst_fill: ", dst_fill);
    dump("dist_func: ", dist_func);
    // Fetches the minimum distance between every source and destination per data.
    isl_pw_qpolynomial *min_dist = minimize_jumps(src_occ, dst_fill, dist_func);
    // First sums cost per dst, then sums cost per dst to get total cost.
    isl_pw_qpolynomial *sum = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(min_dist));
    // Grabs the return value as an isl_val.
    isl_val *sum_extract = isl_pw_qpolynomial_eval(sum, isl_point_zero(isl_pw_qpolynomial_get_domain_space(sum)));
    // Converts isl_val to int.
    long ret = isl_val_get_num_si(sum_extract);

    // Frees val.
    isl_val_free(sum_extract);

    return ret;
}

long analyze_jumps(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func)
{
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
    long ret = analyze_jumps(
        p_src_occupancy,
        p_dst_fill,
        p_dist_func
    );

    // Frees the isl objects.
    isl_ctx_free(p_ctx);

    return ret;
}

/**
 * Analyzes the latency of a memory access by finding the minimum path from
 * every source to every destination for a particular data, then taking the max
 * of all the minimum paths for that data, then taking the max of all the
 * latencies for all the data.
 * 
 * @param __isl_take p_src_occupancy    A map relating source location and the 
 *                                      data occupied.
 * @param __isl_take p_dst_fill         A map relating destination location and
 *                                      the data requested.
 * @param __isl_take dist_func          The distance function to use, as a map.
 */
long analyze_latency (
    isl_map *src_occ, 
    isl_map *dst_fill, 
    isl_map *dist_func
) {
    // Fetches the minimum distance between every source and destination per data.
    isl_pw_qpolynomial *p_min_dist = minimize_jumps(src_occ, dst_fill, dist_func);
    // Computes the maximum of minimum distances for every data.
    isl_val *p_max_min_dist = isl_pw_qpolynomial_max(p_min_dist);
    int ret = isl_val_get_num_si(p_max_min_dist);

    // Frees the isl objects.
    isl_val_free(p_max_min_dist);

    return ret;
}

/**
 * A wrapper for analyze_latency that takes in strings instead of isl objects.
 * 
 * @param src_occupancy     A string representation of a map relating source
 *                          location and the data occupied.
 * @param dst_fill          A string representation of a map relating destination
 *                          location and the data requested.
 * @param dist_func         A string representation of a distance function to use.
 *
 * @return                  A string representation of the maximum latency.
 */
long analyze_latency (
    const std::string& src_occupancy, 
    const std::string& dst_fill, 
    const std::string& dist_func
) {
    // Creates a new isl context.
    isl_ctx *ctx = isl_ctx_alloc();

    // Reads the string representations of the maps into isl objects.
    isl_map *p_src_occ = isl_map_read_from_str(ctx, src_occupancy.c_str());
    isl_map *p_dst_fill = isl_map_read_from_str(ctx, dst_fill.c_str());
    isl_map *p_dist_aff = isl_map_read_from_str(ctx, dist_func.c_str());
    // Calls the isl version of analyze_latency.
    long ret = analyze_latency(p_src_occ, p_dst_fill, p_dist_aff);

    // Frees the isl objects.
    isl_ctx_free(ctx);

    return ret;
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

    // Allocates computer memory for the isl space where dist calculations are done.
    isl_space *p_dist_space = isl_space_alloc(p_ctx, 0, dst_dims.size(), src_dims.size());

    // Programmatically creates and binds the dst and src dimensions to the dist space.
    for (int i = 0; i < src_dims.size(); i++)
    {
        p_dist_space = isl_space_set_dim_id(
            p_dist_space, isl_dim_in, i, isl_id_alloc(p_ctx, dst_dims[i].c_str(), NULL)
        );
        p_dist_space = isl_space_set_dim_id(
            p_dist_space, isl_dim_out, i, isl_id_alloc(p_ctx, src_dims[i].c_str(), NULL)
        );
    }

    // Wraps the space into a set space.
    p_dist_space = isl_space_wrap(p_dist_space);
    // Converts it into a local space.
    isl_local_space *p_dist_local = isl_local_space_from_space(p_dist_space);

    // Total Manhattan distance affine function.
    isl_pw_aff *nd_manhattan_metric = isl_pw_aff_zero_on_domain(
        isl_local_space_copy(p_dist_local)
    );
    // Constructs all the absolute value affines per dimension and adds to metric.
    for (int i = 0; i < dst_dims.size(); i++)
    {
        // Constructs the affine for the src.
        isl_pw_aff *src_aff = isl_pw_aff_var_on_domain(
            isl_local_space_copy(p_dist_local), isl_dim_set, dst_dims.size() + i
        );
        // Constructs the affine for the dst.
        isl_pw_aff *dst_aff = isl_pw_aff_var_on_domain(
            isl_local_space_copy(p_dist_local), isl_dim_set, i
        );

        // Subtracts the dst. 
        isl_pw_aff *p_aff = isl_pw_aff_sub(src_aff, dst_aff);
        // Grabs the negation.
        isl_pw_aff *p_neg_aff = isl_pw_aff_neg(isl_pw_aff_copy(p_aff));
        // Constructs the affine for the absolute value.
        isl_pw_aff *p_abs_aff = isl_pw_aff_max(p_aff, p_neg_aff);

        // Adds the absolute value affine to the vector.
        nd_manhattan_metric = isl_pw_aff_add(nd_manhattan_metric, p_abs_aff);
    }

    // Grabs the return value as a string.
    std::string ret = isl_pw_aff_to_str(nd_manhattan_metric);

    // Frees the isl objects.
    isl_local_space_free(p_dist_local);
    isl_pw_aff_free(nd_manhattan_metric);
    isl_ctx_free(p_ctx);

    return ret;
}

/**
 * Calculates the latency of a memory access on a ring.
 * 
 * @param n     The circumference of the torus. 
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
    isl_space *p_dist_space = isl_space_alloc(p_ctx, 0, 1, 1);
    
    // Sets the dst and src as members of the space.
    p_dist_space = isl_space_set_dim_id(p_dist_space, isl_dim_in, 0, dst_id);
    p_dist_space = isl_space_set_dim_id(p_dist_space, isl_dim_out, 0, src_id);

    // Wraps the space into a set space.
    p_dist_space = isl_space_wrap(p_dist_space);

    // Creates p_dist as a local space.
    isl_local_space *p_dist_local = isl_local_space_from_space(p_dist_space);

    // Creates the src and dst affines.
    isl_pw_aff *src_aff = isl_pw_aff_var_on_domain(isl_local_space_copy(p_dist_local), isl_dim_set, 0);
    isl_pw_aff *dst_aff = isl_pw_aff_var_on_domain(p_dist_local, isl_dim_set, 1);

    // Subtracts the dst from the src.
    isl_pw_aff *src_sub_dst_aff = isl_pw_aff_sub(src_aff, dst_aff);
    // Subtracts the src from the dst.
    isl_pw_aff *dst_sub_src_aff = isl_pw_aff_neg(isl_pw_aff_copy(src_sub_dst_aff));

    // Creates an isl_val for the torus circumference.
    isl_val *p_circumference = isl_val_int_from_si(p_ctx, n);

    // Creates the moduli for both differences.
    isl_pw_aff *src_sub_dst_mod_n_aff = isl_pw_aff_mod_val(src_sub_dst_aff, isl_val_copy(p_circumference));
    isl_pw_aff *dst_sub_src_mod_n_aff = isl_pw_aff_mod_val(dst_sub_src_aff, p_circumference);

    // Combines the moduli into a single piecewise affine.
    isl_pw_aff *p_dist = isl_pw_aff_min(src_sub_dst_mod_n_aff, dst_sub_src_mod_n_aff);
    // Grabs the return value as a string.
    std::string ret = isl_pw_aff_to_str(p_dist);

    // Frees the isl objects.
    isl_pw_aff_free(p_dist);
    isl_ctx_free(p_ctx);

    return ret;
}