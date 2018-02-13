//
// Created by wz on 18-1-16.
//

#ifndef EVALUATOR1_EVALUATOR_H
#define EVALUATOR1_EVALUATOR_H

#include <igraph.h>

#include <queue>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

using namespace std;

void init(igraph_t* g);

void evaluate(string & query_file, string & result_file);

#endif //EVALUATOR1_EVALUATOR_H
