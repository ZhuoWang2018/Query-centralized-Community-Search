//
// Created by wz on 17-9-10.
//
#include <igraph.h>
int get_simple_graph(const char * name, igraph_t * g) {
    FILE *ifile = fopen(name, "r");
    if (ifile == 0) {
        printf("file read failed!\n");
        return 1;
    }
    int errno = igraph_read_graph_edgelist(g, ifile, 0, IGRAPH_UNDIRECTED);
    fclose(ifile);
    if (errno == IGRAPH_PARSEERROR) {
        printf("%s is not a correct graph format~\n", name);
        return 2;
    }
    return 0;
}