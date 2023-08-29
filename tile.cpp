#include "tile.hpp"

int main(int argc, char const *argv[])
{
    /** NECESSARY FOR THE LANGUAGE'S GLOBAL ENVIRONMENT **/
    /** ASSUMED GIVENS FROM USER **/
    // Creates an ISL context.
    isl_ctx *ctx = isl_ctx_alloc();
    /** Creates the topology. **/
    // The SRCs mapping to some unknown occupancy.
    isl_map *src_occ = isl_map_read_from_str(ctx, 
        R"SRCS({ [xs, ys] -> [data] |
            (0 <= xs < 2) and
            (0 <= ys < 2) and
            0 <= data < 16
        })SRCS"
    );
    // The DSTs mapping to a known quantity.
    isl_map *dst_fill = isl_map_read_from_str(ctx, 
        R"DSTS({ [xd, yd] -> [data] |
            (0 <= xd < 4) and
            (0 <= yd < 4) and
            (4yd <= data < 4yd + 4) and
            0 <= data < 16
        })DSTS"
    );

    /** PROGRAMMATIC GENERATION WITH TILE **/
    isl_map *tiling = tile(
        0,
        isl_map_get_space(src_occ),
        8,
        1
    );
    isl_map *subtiling = tile(
        0,
        isl_map_get_space(tiling),
        4,
        0
    );
    isl_map_dump(subtiling);
    src_occ = isl_map_intersect(isl_map_intersect(src_occ, tiling), subtiling);
    isl_map_dump(src_occ);
    isl_map_free(src_occ);
    isl_map_free(dst_fill);

    isl_ctx_free(ctx);

    return 0;
}

/**
 * Creates an ISL set that restricts the data domain to a tiling split along a
 * certain axis.
 * 
 * Read as: Tile the data axis in position data_dim from 
 * src_space in blocks of n consecutive elements along src axis axis_dim.
 * 
 * @param data_dim  __isl_keep  The data axis index.
 * @param src_space __isl_take  The space in which the axis is defined.
 * @param n         __isl_keep  The number of elements in a block.
 * @param axis_dim  __isl_keep  The axis index to tile along.
 */
isl_map *tile(
    int data_dim,
    isl_space *src_space,
    int n, 
    int axis_dim
) {
    /* Creates the tiling restriction */
    // Allocates local space for the tiling restriction.
    isl_local_space *tile_local_space = isl_local_space_from_space(src_space);
    // Creates n*axis <= data (equiv. to data - n*axis >= 0)
    isl_constraint *tile_lower = isl_constraint_alloc_inequality(isl_local_space_copy(tile_local_space));
    tile_lower = isl_constraint_set_coefficient_si(tile_lower, isl_dim_in, axis_dim, -n);
    tile_lower = isl_constraint_set_coefficient_si(tile_lower, isl_dim_out, data_dim, 1);
    // Creates n*axis + n > data (equiv. to n*axis + n - data - 1 >= 0)
    isl_constraint *tile_upper = isl_constraint_alloc_inequality(tile_local_space);
    tile_upper = isl_constraint_set_coefficient_si(tile_upper, isl_dim_in, axis_dim, n);
    tile_upper = isl_constraint_set_coefficient_si(tile_upper, isl_dim_out, data_dim, -1);
    tile_upper = isl_constraint_set_constant_si(tile_upper, n - 1);
    // Creates the tiling restriction.
    isl_basic_map *tile = isl_basic_map_from_constraint(tile_lower);
    tile = isl_basic_map_add_constraint(tile, tile_upper);

    // Returns the tiling restriction
    return isl_map_from_basic_map(tile);
}

/**
 * Creates an ISL set that expands the data data domain in a certain src to duplicate
 * it over a certain axis.
 * 
 * Read as: Replicate the feature n times along the given src axis axis_dim.
 * 
 * @param feature  __isl_take   The binding feature to replicate.
 * @param n        __isl_keep   The number of times to replicate the feature.
 * @param axis_dim __isl_keep   The axis index to replicate along (assumed to be
 *                              an src axis index in the space of feature).
 */
isl_map *replicate(
    isl_map *feature,
    int n,
    int axis_dim
) {
    return nullptr;
}