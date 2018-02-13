//
// Created by wz on 18-1-18.
//

#include "mdc_test.h"

void split_mdc(string str, igraph_vector_int_t * result) {
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


void mdc_test(string & file_name, string & query_file, string & result_file) {
    igraph_t g;
    get_simple_graph(file_name.data(), &g);

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::out);

    int community_counter = 0;

    string s;
    igraph_vector_int_t query_points;
    igraph_vector_int_init(&query_points, 0);
    unordered_set<int> fetched_points;

    chrono::duration<double> elapsed;
    chrono::time_point<chrono::system_clock> start, end;
    init_mdc_ori(&g);
    while(getline(query_f, s)) {
        split_mdc(s, &query_points);
        start = chrono::system_clock::now();;
        fetched_points = find_mdc(&g, &query_points, 100);
        end = chrono::system_clock::now();
        printf("one query has been answered!\n");
        elapsed += (end - start);
        stringstream s_stream;
        auto it_point = fetched_points.begin();
        while(it_point != fetched_points.end()) {
            s_stream<<*it_point<<" ";
            it_point ++;
        }
        string result = s_stream.str();
        result_f<<result.substr(0, result.size() - 1)<<"\n";

        community_counter++;
        if(community_counter % 10 == 0) {
            printf("i am still working, on %d\n", community_counter);
            printf("time elapsed: %fs\n", elapsed.count());
        }
    }
    destroy_mdc_ori();
    igraph_vector_int_destroy(&query_points);
    igraph_destroy(&g);
    query_f.close();
    result_f.close();

}