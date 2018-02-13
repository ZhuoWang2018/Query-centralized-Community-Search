//
// Created by wz on 18-1-30.
//

#ifndef WEIGHTED_DENSITY_EVALUATOR_EVALUATOR_H
#define WEIGHTED_DENSITY_EVALUATOR_EVALUATOR_H

#include <igraph.h>

#include <queue>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

using namespace std;

void init(igraph_t* g);

void evaluate_weighted_density(string & query_file, string & result_file, bool single_included);

#endif //WEIGHTED_DENSITY_EVALUATOR_EVALUATOR_H
