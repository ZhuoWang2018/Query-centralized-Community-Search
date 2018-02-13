//
// Created by wz on 18-1-17.
//
#include "lctc_test.h"

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

void lctc_test(string & graph_file,
               string & query_file,
               string & result_file) {
    igraph_t g;

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::out);

    get_ml_graph(graph_file.data(), &g);

    init_lctc(&g);

    string s;
    int counter = 0;
    int upper_bound = 1000;
    igraph_vector_int_t query_vertices;
    igraph_vector_int_init(&query_vertices, 0);

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed = start - start;

    while(getline(query_f, s)) {
        split_query(s, &query_vertices);
        start = chrono::system_clock::now();


        auto community = fetch_community(&query_vertices, upper_bound);
        end = chrono::system_clock::now();
        stringstream s_stream;
        for(const auto & element: community)
            s_stream<< element <<" ";

        string result = s_stream.str();
        result_f<<result.substr(0, result.size() - 1)<<"\n";

        printf("one query have answered!\n");
        elapsed += (end - start);

        counter++;
        if(counter % 10 == 0) {
            printf("loc-kd is still working, on %d\n", counter);
            printf("time elapsed:%fs\n", elapsed.count());
        }
    }

    igraph_vector_int_destroy(&query_vertices);
    igraph_destroy(&g);
    query_f.close();
    result_f.close();

}
