#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <isl/id.h>
#include <isl/map.h>
#include <isl/set.h>
#include <isl/space.h>

isl_set *tile(
    isl_ctx *ctx, 
    isl_set *data,
    int n, 
    isl_id *axis,
    isl_set *src_space
);