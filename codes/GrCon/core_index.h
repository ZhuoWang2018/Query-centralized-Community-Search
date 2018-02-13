//
// Created by wz on 17-9-10.
//

#ifndef GRCON_CORE_INDEX_H
#define GRCON_CORE_INDEX_H

#include <igraph.h>

#include <queue>
#include <sstream>
#include <fstream>

#include <unordered_map>
#include <unordered_set>

using namespace std;


void core_decomposition(igraph_t * graph, fstream & Hf);

#endif //GRCON_CORE_INDEX_H
