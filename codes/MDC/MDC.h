//
// Created by wz on 18-1-17.
//

#ifndef MDC_MDC_H
#define MDC_MDC_H

#include <igraph.h>

#include <queue>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>

using namespace std;

void init_mdc_ori(igraph_t * g_data);

void destroy_mdc_ori();

unordered_set<int> find_mdc(igraph_t * ori_graph,
                            igraph_vector_int_t * query_vertices, int upper_bound);


#endif //MDC_MDC_H
