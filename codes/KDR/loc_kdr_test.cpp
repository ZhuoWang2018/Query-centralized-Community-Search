//
// Created by wz on 18-2-7.
//

#include "loc_kdr_test.h"

void loc_kdr_test(string & graph_file,
                 string & query_file,
                 string & result_file,
                 int steiner_type, double c, double alpha) {

    igraph_t g;

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::out);

    get_simple_graph(graph_file.data(), &g);

    init_loc(&g, steiner_type);
    set_c(c);
    set_alpha(alpha);

    string s;
    int upper_bound = 200;
    int community_counter = 0;
    igraph_vector_int_t query_vertices;
    igraph_vector_int_init(&query_vertices, 0);

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed = start - start;
    while(getline(query_f, s)) {
        split_query(s, &query_vertices);
        start = chrono::system_clock::now();

        if(igraph_vector_int_size(&query_vertices) == 1) {
            auto communities = get_loc_kdr_single(&query_vertices, upper_bound);
            end = chrono::system_clock::now();
            result_f<<"#\n";
            for(const auto & community: communities) {
                stringstream s_stream;
                for(const auto & element: community)
                    s_stream<<element<<" ";
                string result = s_stream.str();
                result_f<<result.substr(0, result.size() - 1)<<"\n";
            }
            result_f<<"#\n";
        } else {
            auto community = get_loc_kdr(&query_vertices, upper_bound);
            end = chrono::system_clock::now();
            stringstream s_stream;
            for(const auto & element: community)
                s_stream<< element <<" ";

            string result = s_stream.str();
            result_f<<result.substr(0, result.size() - 1)<<"\n";
        }

        printf("one query have answered!\n");
        elapsed += (end - start);

        community_counter++;
        if(community_counter % 10 == 0) {
            printf("loc-kd is still working, on %d\n", community_counter);
            printf("time elapsed:%fs\n", elapsed.count());
        }
    }

    igraph_vector_int_destroy(&query_vertices);
    igraph_destroy(&g);
    query_f.close();
    result_f.close();
}