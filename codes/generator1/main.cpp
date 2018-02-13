#include <iostream>

#include "generator.h"

int main(int argc, char** argv) {

    if(argc < 3) {
        printf("unmatched parameters!\n");
        printf("graph file!\n");
	printf("query file!\n");
        return 0;
    }

    string graph_file(argv[1]);
    string query_dir(argv[2]);

    generate_queries(graph_file, query_dir);

}
