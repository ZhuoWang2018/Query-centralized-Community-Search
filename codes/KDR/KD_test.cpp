//
// Created by wz on 18-2-7.
//

#include "KD_test.h"



void split_query(string str, igraph_vector_int_t * result) {
    string::size_type pos;
    igraph_vector_int_clear(result);
    string pattern = " ";
    str += pattern;
    unsigned long size = (int)str.size();

    for(unsigned long i = 0; i < size; i++) {
        pos = str.find(" ", i);
        if(pos < size) {
            string s = str.substr(i, pos - i);
            int number;
            sscanf(s.data(), "%d", &number);
            igraph_vector_int_push_back(result, number);
            i = pos + pattern.size() - 1;
        }
    }
}

void kd_test(string & graph_file,
                   string & query_file, string & result_file) {

    igraph_t g;

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::out);

    get_ml_graph(graph_file.data(), &g);

    init_kd(&g);

    string s;
    int community_counter = 0;
    unordered_set<int> community;
    igraph_vector_int_t query_vertices;
    igraph_vector_int_init(&query_vertices, 0);

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed = start - start;
    while(getline(query_f, s)) {
        split_query(s, &query_vertices);
        start = chrono::system_clock::now();
        community = get_kd(&query_vertices);
        end = chrono::system_clock::now();
        printf("one query have answered!\n");
        elapsed += (end - start);

        stringstream s_stream;
        for(const auto & element: community)
            s_stream<< element <<" ";

        string result = s_stream.str();
        result_f<<result.substr(0, result.size() - 1)<<"\n";

        community_counter++;
        if(community_counter % 10 == 0) {
            printf("acc plus is still working, on %d\n", community_counter);
            printf("time elapsed:%fs\n", elapsed.count());
        }
    }

    igraph_vector_int_destroy(&query_vertices);
    igraph_destroy(&g);
    query_f.close();
    result_f.close();
}


void kd_offline(string & graph_file, string & graph_out_file) {
    igraph_t g;
    igraph_i_set_attribute_table(&igraph_cattribute_table);
    get_simple_graph(graph_file.data(), &g);
    k_core_decomposition(&g);
    write_ml_graph(graph_out_file.data(), g);
    igraph_destroy(&g);
    printf("k-core decomposition completed!\n");
}