#include <isl/map.h>

int main(int argc, char* argv[])
{
    auto p_ctx = isl_ctx_alloc();
    auto p_src_occupancy = isl_map_read_from_str(
        p_ctx,
        "{ [xs, ys] -> [d0, d1] : d0=x and 0 <= d1 < 8 and 0 <= x < 8 and 0 <= y < 8 }"
    );
    auto p_dst_fill = isl_map_read_from_str(
        p_ctx,
        "{ [xs, ys] -> [d0, d1] : d0=x and d1=y and 0 <= x < 8 and 0 <= y < 8 }"
    );
  
    analyze_latency(p_src_occupancy, p_dst_fill);
}