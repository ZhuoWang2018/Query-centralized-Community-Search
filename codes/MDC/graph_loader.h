//
// Created by wz on 18-1-17.
//

#ifndef MDC_GRAPH_LOADER_H
#define MDC_GRAPH_LOADER_H

#include <igraph.h>

//used for loading and writing graphs to disks

// read simple graph from disks
int get_simple_graph(const char * name, igraph_t * g);

// write graph with attributes into disks "XX.GraphML"
void write_ml_graph(const char * name, igraph_t & g);

// read graph with attributes from disks "XX.GraphML"
int get_ml_graph(const char * name, igraph_t * g);

#endif //MDC_GRAPH_LOADER_H
