#include <iostream>
#include "mdc_test.h"

int main(int argc, char** argv) {
    if(argc < 4) {
        printf("unmatched parameters!\n");
        printf("absolute path for graph file\n");
        printf("absolute path for query file\n");
        printf("absolute path for result file\n");
        return 0;
    }

    string graph_file(argv[1]);
    string query_file(argv[2]);
    string result_file(argv[3]);

    mdc_test(graph_file, query_file, result_file);

    return 0;
}