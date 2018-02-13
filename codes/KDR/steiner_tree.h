//
// Created by wz on 18-2-7.
//

#ifndef KDR_STEINER_TREE_H
#define KDR_STEINER_TREE_H


#include <igraph.h>

#include <cfloat>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/heap/d_ary_heap.hpp>


#include "k_core_decomp.h"

using namespace std;


struct edge{
    int from;
    int to;
    edge(int data_from, int data_to) {
        from = data_from;
        to = data_to;
    }
};

typedef edge* EDGE;

struct edge_hash{
    size_t operator()(const EDGE & e) const {
        return std::hash<int>()(e->from) ^ (std::hash<int>()(e->to) << 1);
    }
};

struct edge_equal{
    bool operator()(const EDGE & le, const EDGE & re) const {
        return le->from == re->from && le->to == re->to;
    }
};

const int arity_number = 4;

struct route{
    double distance;
    int current_vertex_idx;
    int previous_vertex_idx;
    int origin_vertex_idx;
    route(double data_distance, int data_previous_vertex_idx,
          int data_current_vertex_idx, int data_origin_vertex_idx) {
        distance = data_distance;
        current_vertex_idx = data_current_vertex_idx;
        previous_vertex_idx = data_previous_vertex_idx;
        origin_vertex_idx = data_origin_vertex_idx;
    }
};

struct route_cmp{
    bool operator()(const route* a, const route* b) const{
        return a->distance > b->distance;
    }
};

typedef boost::heap::d_ary_heap<route*,
        boost::heap::mutable_<true>,
        boost::heap::arity<arity_number>,
        boost::heap::stable<true>,
        boost::heap::compare<route_cmp>>::handle_type handle_t;

struct pres_node{
    int vertex_idx;
    double prestige;
    double neighbors;
    pres_node(int data_vertex_idx, double data_prestige, double data_neighbors) {
        vertex_idx = data_vertex_idx;
        prestige = data_prestige;
        neighbors = data_neighbors;
    }
};

struct pres_cmp{//the top elements are of bigger prestige
    bool operator()(const pres_node* a, const pres_node* b) const{
        return a->prestige / a->neighbors < b->prestige / b->neighbors;
    }
};

typedef boost::heap::d_ary_heap<pres_node*,
        boost::heap::mutable_<true>,
        boost::heap::arity<arity_number>,
        boost::heap::stable<true>,
        boost::heap::compare<pres_cmp>>::handle_type pres_handle_t;

struct peak_node{
    int prestige_rank = -1;
    int distance_rank = -1;
    double distance;
    double prestige;
};

void init_steiner(igraph_t* g, int type);

void set_steiner_c(double data);

void set_steiner_alpha(double data);

unordered_set<int> get_steiner_tree(igraph_vector_int_t * query_vertices);

#endif //KDR_STEINER_TREE_H
