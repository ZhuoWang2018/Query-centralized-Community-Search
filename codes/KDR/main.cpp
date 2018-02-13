#include <iostream>

#include "KD_test.h"
#include "loc_kdr_test.h"

int main(int argc, char** argv) {
    if(argc < 5) {
        printf("unmatched parameters!\n");
        printf("Algorithms: 1-> KD; 2-> Loc-kdr\n");
        printf("absolute path for graph file\n");
        printf("absolute path for query file\n");
        printf("absolute path for result file\n");
        printf("optional steiner_type: 1-> mehlhorn; 2-> prune; 3-> heuristic(default)\n");
        printf("optional steiner c: default 3.0\n");
        printf("optional steiner alpha: default 0.5\n");
        printf("\n");
        return 0;
    }

    int algorithm_type;
    sscanf(argv[1], "%d", &algorithm_type);

    string graph_file(argv[2]);
    string query_file(argv[3]);
    string result_file(argv[4]);

    int steiner_type = 3;
    double c = 3.0;
    double alpha = 0.5;

    if(argc >= 6)
        sscanf(argv[5], "%d", &steiner_type);

    if(argc >= 7)
        sscanf(argv[6], "%lf", &c);

    if(argc >= 8)
        sscanf(argv[7], "%lf", &alpha);

    if(algorithm_type == 1)
        kd_test(graph_file, query_file, result_file);
    else if(algorithm_type == 2)
        loc_kdr_test(graph_file, query_file, result_file, steiner_type, c, alpha);

    return 0;
}