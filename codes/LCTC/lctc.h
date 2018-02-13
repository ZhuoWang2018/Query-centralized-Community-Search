//
// Created by wz on 18-1-17.
//

#ifndef LCTC_LCTC_H
#define LCTC_LCTC_H

#include "steiner_tree.h"


void init_lctc(igraph_t * g_data);

unordered_set<int> fetch_community(
        igraph_vector_int_t * query_vertices, int upper_bound);
#endif //LCTC_LCTC_H
