#include <iostream>

#include "evaluator.h"
#include "graph_loader.h"

int main(int argc, char** argv) {
    if(argc < 4) {
        printf("unmatched parameters!\n");
        printf("graph file\n");
        printf("query file\n");
        printf("result file\n");
        printf("single included\n");
        return 0;
    }

    string graph_file(argv[1]);
    string query_file(argv[2]);
    string result_file(argv[3]);
    string bool_tag(argv[4]);

    bool single_included = true;
    if(bool_tag == "false")
        single_included = false;


//    string graph_file = "/home/wz/test-0828/amazon/com-amazon.ungraph";
//    string query_file = "/home/wz/test-0828/multiple.query";
//    string result_file = "/home/wz/test-0828/multiple.lctc";
//    bool single_included = false;

    igraph_t graph;
    get_simple_graph(graph_file.data(), &graph);
    init(&graph);
    evaluate_weighted_density(query_file, result_file, single_included);

    igraph_destroy(&graph);

    return 0;
}