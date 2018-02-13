//
// Created by wz on 18-1-16.
//

#ifndef GENERATOR1_GENERATOR_H
#define GENERATOR1_GENERATOR_H

#include <queue>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>

#include "graph_loader.h"

using namespace std;
using namespace boost;

void generate_queries(string & file_name,
                      string & query_directory);

#endif //GENERATOR1_GENERATOR_H
