//
// Created by wz on 18-2-7.
//

#ifndef KDR_KD_TEST_H
#define KDR_KD_TEST_H

#include <chrono>
#include <fstream>
#include <sstream>

#include "KD.h"
#include "graph_loader.h"


void split_query(string str, igraph_vector_int_t * result);

void kd_test(string & graph_file,
                   string & query_file, string & result_file);


#endif //KDR_KD_TEST_H
