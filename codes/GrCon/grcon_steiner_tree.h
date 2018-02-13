//
// Created by wz on 17-9-12.
//

#ifndef GRCON_GRCON_STEINER_TREE_H
#define GRCON_GRCON_STEINER_TREE_H

#include <queue>
#include <igraph.h>
#include <unordered_map>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>
#include <boost/heap/d_ary_heap.hpp>

using namespace std;

const int arity_number = 4;

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

struct route{
    int min_truss;
    double distance;
    int current_vertex_idx;
    int previous_vertex_idx;
    int origin_vertex_idx;
    route(int data_min_truss, double data_distance, int data_previous_vertex_idx,
          int data_current_vertex_idx, int data_origin_vertex_idx) {
        min_truss = data_min_truss;
        distance = data_distance;
        current_vertex_idx = data_current_vertex_idx;
        previous_vertex_idx = data_previous_vertex_idx;
        origin_vertex_idx = data_origin_vertex_idx;
    }
};

struct edge_equal{
    bool operator()(const EDGE & le, const EDGE & re) const {
        return le->from == re->from && le->to == re->to;
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

//void initialize_steiner_tree(igraph_t* gtemp);

void initialize_steiner_tree(igraph_t* gtemp,
                             boost::dynamic_bitset<> * Htemp);


unordered_set<int> get_steiner_tree(igraph_vector_int_t * vertex_idx_array);

#endif //GRCON_GRCON_STEINER_TREE_H
