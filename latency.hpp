#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Includes ISL affine list/piecewise functions.
#include <isl/aff.h>
#include <isl/ilp.h>
// Includes isl qpolynomials.
#include <isl/polynomial.h>
// Includes ISL maps/binary relations.
#include <isl/map.h>
// Includes ISL ids and dspaces.
#include <isl/id.h>
#include <isl/space.h>
// Imports ISL sets.
#include <isl/set.h>
// Imports ISL val.
#include <isl/val.h>
#include <barvinok/isl.h>
#include <barvinok/polylib.h>

long analyze_jumps(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_pw_aff *dist_func);
long analyze_jumps(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);
long analyze_latency(isl_map *p_src_occupancy, isl_map *p_dst_fill, isl_pw_aff *dist_func);
long analyze_latency(const std::string& src_occupancy, const std::string& dst_fill, const std::string& dist_func);
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims);
std::string n_long_ring_metric(long n);

// Defines debug variables from environment variables.
#include <string.h>
bool islIntermediates = (getenv("ISL_INTERMEDIATES") != NULL) &&
                        (strcmp(getenv("ISL_INTERMEDIATES"), "0") != 0);

// Defines a function to programatically generate an n-dimensional Manhattan distance function.
std::string nd_manhattan_metric(std::vector<std::string> src_dims, std::vector<std::string> dst_dims);
// Defines a function to programatically generate an n-circumference ring distance function.
std::string n_long_ring_metric(long ring_circumference);

// Defines debug dump function.
void dump(std::string str, isl_map *map)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_map_dump(map);
    }
}

void dump(std::string str, isl_pw_aff *pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_pw_aff_dump(pw_aff);
    }
}

void dump(std::string str, isl_multi_pw_aff *multi_pw_aff)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_pw_aff_dump(multi_pw_aff);
    }
}

void dump(std::string str, isl_multi_val *multi_val)
{
    if (islIntermediates)
    {
        std::cout << str << std::endl;
        isl_multi_val_dump(multi_val);
    }
}

#define DUMP(varname) dump(#varname, varname)

/// @brief Strings representing the src and dst datum holds/requests in ISL.
struct binding_struct
{
    const std::string srcs;
    const std::string dsts;
};
typedef std::unique_ptr<binding_struct> binding;
/** 
 * @brief Defines the struct that comprises the result of folding and the unique
 * ptr to it that represents what is returned by fold.
 */
struct fold_struct
{
    const long cost;
    const std::string folded_repr;
};
typedef std::unique_ptr<fold_struct> fold_result;
/// @brief Defines the struct characterizing the collapsing behavior of a layer.
struct collapse_struct
{
    const std::string src_collapser;
    const std::string dst_collapser;
};
typedef std::shared_ptr<collapse_struct> collapse;

/// @brief Virtual class that represents a computation Layer for costs.
class Layer {
    public:
        /// @brief The cost formula of the folding step for this layer.
        const std::string crease_costs;
        /** 
         * @brief The folding action to a multicastable representation after
         * calculating the cost of folding. 
         */
        const std::string fold_formula;
        /// @brief The cost formula of the multicasting step for this layer.
        const std::string multicast_costs;
        /// @brief The src collapse formulation for the next layer.
        const std::string src_collapser;
        /// @brief The dst collapse formulation for the next layer.
        const std::string dst_collapser;
        /// @brief The context the layer is in.
        isl_ctx *const ctx;
    // public:
    //     virtual binding evaluate
};

// class Mesh: public Layer {

// }