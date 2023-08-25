#include <string>
#include <unordered_map>
#include <vector>

#include <isl/id.h>
#include <isl/map.h>
#include <isl/set.h>
#include <isl/space.h>

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

    // Dumps the maps.
    isl_map_dump(src_occ);
    isl_map_dump(dst_fill);
    isl_set_dump(data_domain);
    /** PROGRAMMATIC GENERATION WITH TILE **/
    // Tile("d", 2, "ys")
    

    return 0;
}

