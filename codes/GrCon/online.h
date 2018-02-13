//
// Created by wz on 17-9-10.
//

#ifndef GRCON_ONLINE_H
#define GRCON_ONLINE_H

#include <igraph.h>

#include <chrono>
#include <vector>
#include <fstream>
#include <sstream>

#include <unordered_set>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/heap/d_ary_heap.hpp>

#include "grcon_steiner_tree.h"
using namespace std;



//struct vital_info {
//    int alpha;
//    int p_t;
//    int p_tt_positive;
//    int p_tt_negative;
//};
//
//struct pscore {
//    int vertex;
//    int connection_score;
//    int degree_score;
//
//    int pp_positive;
//    int pp_negative;
//    pscore(int data_vertex, int data_connection_sore, int data_degree_score) {
//        vertex = data_vertex;
//        connection_score = data_connection_sore;
//        degree_score = data_degree_score;
//    }
//};
//
//struct pscore_cmp {
//    bool operator()(const pscore* a, const pscore* b) const {
//        if(a->connection_score < b->connection_score)
//            return true;
//        else if(a->connection_score == b->connection_score
//           && a->degree_score < b->degree_score)
//            return true;
//        return false;
//    }
//};
//
//struct deg {
//    int degree;
//    deg(int data_degree) {
//        degree = data_degree;
//    }
//};
//
//struct deg_cmp {
//    bool operator()(const deg* a, const deg* b) const {
//        return a->degree > b->degree;
//    }
//};
//
//typedef boost::heap::d_ary_heap<pscore*,
//        boost::heap::mutable_<true>,
//        boost::heap::arity<arity_number>,
//        boost::heap::stable<true>,
//        boost::heap::compare<pscore_cmp>>::handle_type handle_score_t;
//
//typedef boost::heap::d_ary_heap<deg*,
//        boost::heap::mutable_<true>,
//        boost::heap::arity<arity_number>,
//        boost::heap::stable<true>,
//        boost::heap::compare<deg_cmp>>::handle_type handle_degree_t;

void grcon_test(string & file_name,
                               string & index_name,
                               string & query_file,
                               string & result_file);


struct NodeRank {
    int node;
    int solutionDegree;

    int rank1;
    int rank2;

    NodeRank(int t_node, int t_min_core,
             unordered_map<int, int> & t_node_scores,
             unordered_set<int> & explicitNeighbors) {
        node = t_node;

        unordered_set<int> neighbors;
        for(const auto & element: explicitNeighbors)
            if(t_node_scores.find(element)
                    != t_node_scores.end())
                neighbors.insert(element);

        solutionDegree = neighbors.size();

        rank1 = 0;
        for(const auto & neighbor: neighbors) {
            if(t_node_scores[neighbor] > 0)
                rank1++;
        }

        rank2 = t_min_core - solutionDegree;
        if(rank2 < 0)
            rank2 = 0;
    }

    void updateRankForSaturation() {
        rank1--;
    }

    void updateRanksForAddition(int addedNode, int t_min_core, unordered_map<int,int> & t_node_scores) {
        if(t_node_scores[addedNode] > 0)
            rank1++;

        solutionDegree++;

        rank2 = t_min_core - solutionDegree;
        if(rank2 < 0)
            rank2 = 0;
    }

    int getRank() {
        return rank1 - rank2;
    }
};

struct NodeRank_cmp {
    bool operator()(const NodeRank & a,  const NodeRank & b) const {
        return (a.rank1 - a.rank2) <= (b.rank1 - b.rank2);
    }
};


typedef boost::heap::d_ary_heap<NodeRank,
        boost::heap::mutable_<true>,
        boost::heap::arity<arity_number>,
        boost::heap::stable<true>,
        boost::heap::compare<NodeRank_cmp>>::handle_type handle_score_t;
#endif //GRCON_ONLINE_H
