//
// Created by wz on 18-1-17.
//

#ifndef LCTC_LCTC_TEST_H
#define LCTC_LCTC_TEST_H

#include <chrono>
#include <string>
#include <fstream>
#include <sstream>

#include "lctc.h"
#include "graph_loader.h"

using namespace std;

void lctc_test(string & graph_file,
               string & query_file,
               string & result_file);

#endif //LCTC_LCTC_TEST_H
