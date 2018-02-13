#include <iostream>

#include "online.h"
#include "core_index.h"
#include "graph_creator.h"

int main(int argc, char** argv) {

    if(argc < 5) {
        printf("unmatched parameters!\n");
        printf("absolute path for graph file\n");
        printf("absolute path for index file\n");
        printf("absolute path for query file\n");
        printf("absolute path for result file\n");
        return 0;
    }

    string graph_file(argv[1]);
    string index_file(argv[2]);
    string query_file(argv[3]);
    string result_file(argv[4]);

//    string graph_file = "/home/wz/test-0828/amazon/com-amazon.ungraph";
//    string index_file = "/home/wz/PycharmProjects/graphsets/amazon/core_H";
//    string query_file = "/home/wz/test-0828/amazon/multiple.query";
//    string result_file = "/home/wz/test-0828/amazon/multiple.res.grcon";

    grcon_test(graph_file, index_file, query_file, result_file);
    return 0;
}