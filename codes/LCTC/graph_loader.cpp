//
// Created by wz on 18-1-17.
//

#include "graph_loader.h"



int get_ml_graph(const char * name, igraph_t * g) {
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    FILE *ifile = fopen(name, "r");
    if (ifile == 0) {
        printf("file read failed\n");
        return 1;
    }
    int errno = igraph_read_graph_graphml(g, ifile, 0);
    fclose(ifile);
    if (errno == IGRAPH_UNIMPLEMENTED || errno == IGRAPH_PARSEERROR) {
        printf("%s is not a correct graph format~\n", name);
        return 2;
    }
    return 0;
}

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

void write_ml_graph(const char * name, igraph_t & g) {
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    FILE * ofile = fopen(name, "w");
    if(ofile == 0) {
        printf("ml file write failed!\n");
        return;
    }
    int errno = igraph_write_graph_graphml(&g, ofile, 0);
    fclose(ofile);

    if(errno == IGRAPH_EFILE) {
        printf("write ml graph failed!\n");
        return;
    }
}