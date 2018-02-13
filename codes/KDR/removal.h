//
// Created by wz on 18-2-7.
//

#ifndef KDR_REMOVAL_H
#define KDR_REMOVAL_H

#include <igraph.h>

#include <queue>
#include <unordered_set>
#include <unordered_map>

using namespace std;

// input: vertices to delete (L)
// output: the vertices need to remove for maintaining k-core
unordered_set<int> maintain_k_core(igraph_t* g,
                                   unordered_set<int>* vertices_to_del,
                                   unordered_set<int>* removed_vertices,
                                   igraph_vector_t* degrees, int k_core);

// find the vertices of largest query distance
// if bulk = true, return the vertices together
// if bulk = false, return only one of the vertices
unordered_set<int> get_max_query_dist_vertices(igraph_t* g,
                                               unordered_set<int>* query_nodes,
                                               unordered_set<int>* removed_vertices,
                                               int & max_distance, bool bulk);

#endif //KDR_REMOVAL_H
