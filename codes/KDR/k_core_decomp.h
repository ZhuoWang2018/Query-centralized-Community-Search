//
// Created by wz on 18-2-7.
//

#ifndef KDR_K_CORE_DECOMP_H
#define KDR_K_CORE_DECOMP_H
#include <igraph.h>


#include <queue>
#include <vector>
#include <unordered_set>
#include <unordered_map>


using namespace std;

int get_induced_vertex(int vertex_name,
                       unordered_map<int, int> * induced_map);

// obtain a maximal connected k-core subgraph Gk from G with specific k
// if output = false, Gk does not exist
bool compute_maximal_k_core_from_ori(unordered_set<int> & query_nodes,
                                     igraph_vector_int_t* query_vertices,
                                     int k_core, igraph_t* g,
                                     igraph_t* k_core_graph, vector<int>* k_core_names,
                                     unordered_set<int> & query_nodes_in_k_core);

bool check_query_connectivity(igraph_t* g, int start_vertex, unordered_set<int>* query_nodes,
                              unordered_set<int>* removed_vertices, unordered_set<int>* pending_vertices);

unordered_set<int> get_query_component(igraph_t* g, int start_vertex,
                                       unordered_set<int>* removed_vertices);

// k-core decomposition for each k in G(V,E)
void k_core_decomposition(igraph_t* g);

// for acc+, find the maximal connected k-core subgraph containing Q with largest k
void compose_maximal_core(igraph_t* g,
                          int & cur_cnum, int & next_cnum,
                          unordered_set<int> & query_nodes,
                          igraph_vector_int_t* query_vertices,
                          igraph_t* k_core_graph, vector<int>* k_core_names,
                          unordered_set<int> & query_nodes_in_k_core);

// find a connected k-core subgraph containing Q with smaller k
void compose_k_core(igraph_t* g, int cnum, int & next_cnum,
                    igraph_vector_int_t* query_vertices,
                    igraph_t* k_core_graph, vector<int>* k_core_names,
                    unordered_set<int> & query_nodes_in_k_core);

void compute_maximal_core(
        unordered_set<int>* query_nodes_in_cur,
        int* maximal_core, igraph_t* graph, vector<int>* names,
        igraph_t* maximal_core_graph, vector<int>* maximal_core_names,
        unordered_set<int> & query_nodes_in_maximal_core);

unordered_map<int, int> induce_subgraph_from_ori(igraph_t* g, unordered_set<int>* vertices,
                                                 igraph_t* induced_graph, vector<int>* induced_names);
#endif //KDR_K_CORE_DECOMP_H
