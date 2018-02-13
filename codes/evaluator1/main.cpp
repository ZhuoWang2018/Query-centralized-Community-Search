#include <iostream>

#include "evaluator.h"
#include "graph_loader.h"


using namespace std;


int main(int argc, char ** argv) {

    if(argc < 4) {
        printf("unmatched parameters!\n");
        printf("graph file\n");
        printf("query file\n");
        printf("result file\n");
        return 0;
    }

    string graph_file(argv[1]);
    string query_file(argv[2]);
    string result_file(argv[3]);

//    string graph_file = "/home/wz/test-0828/amazon/com-amazon.ungraph";
//    string result_file = "/home/wz/test-0828/vary_number.1.res.loc-kd";

//    string graph_file = "/home/wz/test-0828/sample";
//    string query_file = "/home/wz/test-0828/query";
//    string result_file = "/home/wz/test-0828/sample.res.acc-plus";

    igraph_t graph;
    get_simple_graph(graph_file.data(), &graph);
    init(&graph);
    evaluate(query_file, result_file);
    igraph_destroy(&graph);

    return 0;
}