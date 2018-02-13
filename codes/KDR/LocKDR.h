//
// Created by wz on 18-2-7.
//

#ifndef KDR_LOCKDR_H
#define KDR_LOCKDR_H

#include "steiner_tree.h"
#include "removal.h"

struct relNode {
    int vertex;
    double relevance;
    int distance;
    relNode(int data_vertex, double data_rel, int data_dist) {
        vertex = data_vertex;
        relevance = data_rel;
        distance = data_dist;
    }
};

struct relNode_can_cmp {
    bool operator()(const relNode & a, const relNode & b) const{
        return a.relevance < b.relevance;
    }
};

struct relNode_bi_cmp {
    bool operator()(const relNode & a, const relNode & b) const {
        if(a.distance < b.distance)
            return true;
        if(a.distance > b.distance)
            return false;
        return a.relevance > b.relevance;
    }
};

void init_loc(igraph_t* graph, int steiner_type);

void set_c(double c);

void set_alpha(double alpha);

vector<unordered_set<int>> get_loc_kdr_single(
        igraph_vector_int_t* query_vertices, int upper_bound);

unordered_set<int> get_loc_kdr(
        igraph_vector_int_t* query_vertices, int upper_bound);

#endif //KDR_LOCKDR_H
