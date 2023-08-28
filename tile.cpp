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
            (0 <= xs < 4 and xs % 2 = 0) and
            (0 <= ys < 4 and ys % 2 = 0)
        })SRCS"
    );
    // The DSTs mapping to a known quantity.
    isl_map *dst_fill = isl_map_read_from_str(ctx, 
        R"DSTS({ [xd, yd] -> [data] |
            (0 <= xd < 4 and xd % 2 = 0) and
            (0 <= yd < 4 and yd % 2 = 0) and
            (4yd <= data < 4yd + 4)
        })DSTS"
    );
    // The data space as a set.
    isl_set *data_domain = isl_set_read_from_str(ctx, 
        "{ [data] | 0 <= data < 16 }"
    );

    /** PROGRAMMATIC GENERATION WITH TILE **/
    tile(
        ctx,
        isl_set_copy(data_domain),
        2,
        isl_id_alloc(ctx, "ys", NULL),
        isl_map_domain(isl_map_copy(src_occ))
    );
    

    return 0;
}

/**
 * Creates an ISL set that restricts the data domain to a tiling split along a
 * certain axis.
 * 
 * Read as: In context ctx, tile the data in blocks of n (consecutive) along
 * this given axis in the given src_space.
 * 
 * @param ctx       __isl_keep  The ISL context.
 * @param data      __isl_take  The data domain being tiled.
 * @param n         __isl_keep  The number of blocks (consecutive elements).
 * @param axis      __isl_take  The axis along which to tile.
 * @param src_space __isl_take  The space in which the axis is defined.
 */
isl_set *tile(
    isl_ctx *ctx, 
    isl_set *data,
    int n, 
    isl_id *axis,
    isl_set *src_space
) {
    /* Calculates the block size of the tiling. */
    // Calculates the number of data elements.
    /// @todo Find a constant-time way to do this.
    uint32_t data_len = 0;
    isl_set_foreach_point(data, 
    [](isl_point *point, void *user) -> isl_stat {
        uint32_t *data_len = (uint32_t *)user;
        data_len[0]++;
        return isl_stat_ok;
    }, &data_len);
    std::cout << "data_len: " << data_len << std::endl;
    // Calculates the block size with ceiling division.
    uint32_t block_size = data_len / n + (data_len % n != 0);
    std::cout << "block_size: " << block_size << std::endl;


    // Calculates the block size of the tiling.
    return nullptr;
}
