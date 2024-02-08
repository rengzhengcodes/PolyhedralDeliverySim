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
    // Defines the src occupancy map as a string.
    std::string src_occupancy = "{[xs, ys] -> [a, b] : a=xs and b=ys and 0 <= xs < 8 and 0 <= ys < 8}";
    // Defines the dst fill map as a string.
    std::string dst_fill = "{[xd, yd] -> [a, b] : a=xd and 0 <= xd < 8 and 0 <= yd < 8 and 0 <= a < 8 and 0 <= b < 8}";

    // Defines the torus circumference.
    int torus_circumference = 8;
    // Defines the distance function string.
    std::string dist_func_str = nd_manhattan_metric({"xs", "ys"}, {"xd", "yd"});
    std::cout << dist_func_str << std::endl;
  
    // long latency = analyze_latency(src_occupancy, dst_fill, dist_func_str);
    // std::cout << "latency: " << latency << std::endl;
    long jumps = analyze_jumps(src_occupancy, dst_fill, dist_func_str);
    std::cout << "jumps:\t " << jumps << std::endl;
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
isl_pw_qpolynomial_fold *minimize_jumps(
    isl_map *p_src_occupancy, 
    isl_map *p_dst_fill, 
    isl_pw_aff *dist_func
) {
    std::cout << "p_dst_fill" << std::endl;
    isl_map_dump(p_dst_fill);
    /* Inverts dst_fill such that data implies dst.
     * i.e. {[xd, yd] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *dst_fill_inverted = isl_map_reverse(
        isl_map_copy(p_dst_fill)
    );
    std::cout << "dst_fill_inv" << std::endl;
    isl_map_dump(dst_fill_inverted);
    /* Inverts src_occupancy such that data implies source.
     * i.e. {[xs, ys] -> [d0, d1]} becomes {[d0, d1] -> [xs, ys]} */
    isl_map *src_occupancy_inverted = isl_map_reverse(p_src_occupancy);
    std::cout << "src_occupancy_inverted" << std::endl;
    isl_map_dump(src_occupancy_inverted);
    /* Takes the factored range of src_occupancy_inverted and dst_fill_inverted
     * to get {[d0, d1] -> [[xd, yd] -> [xs, ys]]} */
    isl_map *data_TO_dst_to_src = isl_map_range_product(
            dst_fill_inverted, src_occupancy_inverted
    );
    std::cout << "data_TO_dst_to_src" << std::endl;
    isl_map_dump(data_TO_dst_to_src);
    /* Wraps dst fill such that the binary relation implies data.
     * i.e. {[[xd, yd] -> [d0, d1]] -> [d0, d1]} */
    isl_map *dst_fill_wrapped = isl_map_range_map(
        isl_map_copy(p_dst_fill)
    );
    std::cout << "dst_fill_wrapped" << std::endl;
    isl_map_dump(dst_fill_wrapped);

    /* Composites dst_fill_wrapped and data_to_dst_to_src to get
     * {[[xd, yd] -> [d0, d1]] -> [[xd', yd'] -> [xs, ys]]} */
    isl_map *dst_to_data_TO_dst_to_src = isl_map_apply_range(
        dst_fill_wrapped, data_TO_dst_to_src
    );

    // Restricts the range such that xd' = xd and yd' = yd.
    for (int i = 0; i < isl_map_dim(p_dst_fill, isl_dim_in); i++)
    {
        /* Restricts the ith element of the output by equating it to the ith
         * element of the input. Treats input and output as if it were flat. */
        dst_to_data_TO_dst_to_src = isl_map_equate(
            dst_to_data_TO_dst_to_src,
            isl_dim_in, i,
            isl_dim_out, i
        );
    };
    isl_map_free(p_dst_fill);

    // Converts the distance function into a pw_qpolynomial.
    isl_pw_qpolynomial *dist_func_pw = isl_pw_qpolynomial_from_pw_aff(dist_func);
    std::cout << "dist func pw:" << isl_pw_qpolynomial_to_str(dist_func_pw) << std::endl;
    // Converts the pw_qpolynomial into a pw_qpolynomial_fold.
    isl_pw_qpolynomial_fold *dist_func_fold = isl_pw_qpolynomial_fold_from_pw_qpolynomial(
        isl_fold_max, dist_func_pw
    );
    
    /* Computes the manhattan distance between the destination for a data and
     * a source for that data. */
    isl_bool b = isl_bool_true;
    std::cout << "dst_to_data_TO-dst_to_src:" << isl_map_to_str(dst_to_data_TO_dst_to_src) << std::endl;
    isl_pw_qpolynomial_fold *manhattan_distance = isl_map_apply_pw_qpolynomial_fold(
        dst_to_data_TO_dst_to_src, dist_func_fold, &b
    );
    assert(b == isl_bool_true);
    // std::cout << "isl_pw_qpolynomial_fold:" << isl_pw_qpolynomial_fold_list_to_str(isl_pw_qpolynomial_fold_to_list(manhattan_distance))
    // << std::endl;

    return manhattan_distance;
}

/**
 * Analyzes the total jumps that a process takes given the source, destination,
 * and distance function.
 * 
 * @param __isl_take p_src_occupancy    A map relating source location and the
 *                                     data occupied.
 * @param __isl_take p_dst_fill         A map relating destination location and
 *                                     the data requested.
 * @param __isl_take dist_func          The distance function to use, as a map.
 */
long analyze_jumps(isl_map *p_src_occ, isl_map *p_dst_fill, isl_pw_aff *p_dist_func)
{
    // prints out inputs
    dump("src_occupancy: ", p_src_occ);
    dump("dst_fill: ", p_dst_fill);
    dump("dist_func: ", p_dist_func);

    // Fetches the minimum distance between every source and destination per data.
    isl_pw_qpolynomial_fold *p_min_dist = minimize_jumps(p_src_occ, p_dst_fill, p_dist_func);

    // Goes over all the qpolynomial_folds, minimizes them, and adds them to the total.
    isl_val *p_total_jumps = isl_val_zero(isl_pw_qpolynomial_fold_get_ctx(p_min_dist));

    // Turns 
    isl_pw_qpolynomial *un_fold = gather_pw_qpolynomial_from_fold(p_min_dist);

    // 
    isl_pw_qpolynomial *sum = isl_pw_qpolynomial_sum(isl_pw_qpolynomial_sum(un_fold));

    // Grabs the return value as a int.
    long ret = isl_val_get_num_si(p_total_jumps);

    // Frees the isl objects.
    isl_val_free(p_total_jumps);

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
    isl_pw_aff *p_dist_func = isl_pw_aff_read_from_str(
        p_ctx,
        dist_func.c_str()
    );

    // prints out inputs
    dump("src_occupancy: ", p_src_occupancy);
    dump("dst_fill: ", p_dst_fill);
    dump("dist_func: ", p_dist_func);
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
    isl_map *p_src_occ, 
    isl_map *p_dst_fill, 
    isl_pw_aff *p_dist_func
) {
    // Fetches the minimum distance between every source and destination per data.
    isl_pw_qpolynomial_fold *p_min_dist = minimize_jumps(p_src_occ, p_dst_fill, p_dist_func);
    // Computes the maximum of minimum distances for every data.
    isl_val *p_max_min_dist = isl_pw_qpolynomial_fold_max(p_min_dist);
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
    isl_pw_aff *p_dist_aff = isl_pw_aff_read_from_str(ctx, dist_func.c_str());
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