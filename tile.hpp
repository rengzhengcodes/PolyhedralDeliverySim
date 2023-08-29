#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <isl/aff.h>
#include <isl/constraint.h>
#include <isl/id.h>
#include <isl/map.h>
#include <isl/set.h>
#include <isl/space.h>

isl_map *tile(
    int data_dim,
    isl_space *src_space,
    int n, 
    int axis_dim
);
isl_map *replicate(
    isl_map *feature,
    int n,
    int axis_dim
);