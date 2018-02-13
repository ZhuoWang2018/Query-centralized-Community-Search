//
// Created by wz on 18-2-7.
//

#ifndef KDR_KD_H
#define KDR_KD_H



#include "removal.h"
#include "k_core_decomp.h"

void init_kd(igraph_t* graph);

unordered_set<int> get_kd(igraph_vector_int_t* query_vertices);


#endif //KDR_KD_H
