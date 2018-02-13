//
// Created by wz on 18-1-16.
//

#include "generator.h"

bool contains(vector<unordered_set<int>> & queries,
              unordered_set<int> query) {
    for(int i = 0; i < queries.size(); i++) {
        if(queries[i] == query)
            return true;
    }
    return false;
}

void get_sorted_degrees(igraph_t* graph,
                        igraph_vector_t & degrees,
                        vector<int> & degree_arrays,
                        vector<int> & degree_starts,
                        vector<int> & degree_ends) {
    igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    degree_arrays.resize((unsigned long)igraph_vcount(graph));

    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    degree_starts.resize((unsigned long)(max_degree + 1));
    degree_ends.resize((unsigned long)(max_degree + 1));

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(graph); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        degree_ends[vertex_degree]++;
    }
}

dynamic_bitset<> bfs_hits(igraph_t* graph, int start_bfs_idx,
                          dynamic_bitset<> & candidates, int inner_distance) {

    dynamic_bitset<> res(candidates.size());

    res.set((unsigned long)start_bfs_idx);

    queue<pair<int, int>*> bfs_queue;
    auto pair_element = new pair<int, int>(start_bfs_idx, 0);
    bfs_queue.push(pair_element);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(!bfs_queue.empty()) {
        auto bfs_e = bfs_queue.front();
        bfs_queue.pop();
        int source = bfs_e->first;
        int distance = bfs_e->second;
        delete bfs_e;

        igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
        for(int k = 0; k < igraph_vector_size(&neighbors); k++) {
            auto neighbor = (int)(VECTOR(neighbors)[k]);
            if(!res.test((unsigned long)neighbor) && distance + 1 <= inner_distance) {
                res.set((unsigned long)neighbor);
                pair_element = new pair<int, int>(neighbor, distance + 1);
                bfs_queue.push(pair_element);
            }
        }
    }

    int dd = (int)res.count();
    igraph_vector_destroy(&neighbors);
    return res & candidates;
}

unordered_set<int> vary_q_number(igraph_t* graph, int inner_distance, int query_number, dynamic_bitset<> & candidates) {

    int order = (int)(rand() % candidates.count());
    auto start_vertex_idx = candidates.find_first();
    for(int i = 0; i < order; i++)
        start_vertex_idx = candidates.find_next(start_vertex_idx);

    if(query_number == 1) {
        unordered_set<int> res;
        int vertex = (int)start_vertex_idx;
        res.insert(vertex);
        return res;
    }

    auto base = bfs_hits(graph, (int)start_vertex_idx, candidates, inner_distance);

    order = (int)(rand() % base.count());
    auto second_vertex_idx = base.find_first();
    for(int i = 0; i < order; i++)
        second_vertex_idx = base.find_next(second_vertex_idx);

    if(second_vertex_idx == start_vertex_idx)
        return unordered_set<int>();

    auto second_base = bfs_hits(graph, (int)second_vertex_idx, candidates, inner_distance);

    auto third_base = base & second_base;

    int third_count = (int)third_base.count();
    if(third_count < query_number)
        return unordered_set<int>();

    unordered_set<int> queries;
    int element = (int)start_vertex_idx;
    queries.insert(element);
    element = (int)second_vertex_idx;
    queries.insert(element);
    while(queries.size() < query_number) {
        order = rand() % third_count;
        auto third_vertex_idx = third_base.find_first();
        for(int i = 0; i < order; i++)
            third_vertex_idx = third_base.find_next(third_vertex_idx);
        int e = (int)third_vertex_idx;
        if(queries.find(e) == queries.end())
            queries.insert(e);
    }

    return queries;
}

void generate_queries_vary_degree(igraph_t* graph, string & query_dir) {
    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    vector<int> degree_arrays;
    vector<int> degree_starts;
    vector<int> degree_ends;
    get_sorted_degrees(graph, degrees, degree_arrays, degree_starts, degree_ends);
    int counter = 0;

    for(double percent = 0.2; percent <= 1.0; percent += 0.2) {

        counter++;

        stringstream dir_stram;
        dir_stram<<query_dir<<"/vary_degree.query."<<(100 - 20*counter);
        fstream query_f(dir_stram.str(), ios::out);

        dynamic_bitset<> candidates((unsigned long)igraph_vcount(graph));
        int start_index = (igraph_vector_size(&degrees) - degree_starts[1]) * (percent - 0.2) + degree_starts[1];
        int end_index = (igraph_vector_size(&degrees) - degree_starts[1]) * percent + degree_starts[1];

        for(int i = start_index; i < end_index; i++)
            candidates.set((unsigned long)degree_arrays[i]);

        printf("candidates size:%d\n", (int)candidates.count());

        vector<unordered_set<int>> queries;
        while(queries.size() < 100) {
            unordered_set<int> query = vary_q_number(graph, 2, 3, candidates);
            if(query.size() > 0 && !contains(queries, query)) {
                queries.push_back(query);
                printf("query size:%d\n", (int)queries.size());
            }
        }
        printf("bucket %f, enough~\n", percent);

        for(int k = 0; k < queries.size(); k++) {

            stringstream s_stream;
            unordered_set<int> sample = queries[k];
            for(const auto & element: sample) {
                printf("degree: %d\n", (int)VECTOR(degrees)[element]);
                s_stream << element << " ";
            }
            string s = s_stream.str();
            query_f<<s.substr(0, s.size() - 1)<<"\n";

        }

        query_f.close();
    }
    igraph_vector_destroy(&degrees);
}

void generate_queries_vary_qnumber(igraph_t* graph, string & query_dir) {

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    vector<int> degree_arrays;
    vector<int> degree_starts;
    vector<int> degree_ends;
    get_sorted_degrees(graph, degrees, degree_arrays, degree_starts, degree_ends);

    dynamic_bitset<> candidates((unsigned long)igraph_vcount(graph));

    int index = (igraph_vector_size(&degrees) - degree_starts[1]) * 0.2 + degree_starts[1];

    int degree_start = degree_starts[VECTOR(degrees)[degree_arrays[index]]];
    int degree_end = degree_ends[VECTOR(degrees)[degree_arrays[index]]];

    for(int i = degree_start; i <= degree_end; i++) {
        candidates.set((unsigned long) (degree_arrays[i]));
    }

    printf("candidates size:%d\n", (int)candidates.count());

    int counter = 0;

    for(int i = 1; i <= 16; i*=2) {

        counter++;

        stringstream dir_stram;
        dir_stram<<query_dir<<"/vary_number.query."<<i;

        fstream query_f(dir_stram.str(), ios::out);

        vector<unordered_set<int>> queries;
        while(queries.size() < 100) {
            unordered_set<int> query = vary_q_number(graph, 2, i, candidates);
            if(query.size() > 0 && !contains(queries, query)) {
                queries.push_back(query);
                printf("query size:%d\n", (int)queries.size());
            }
        }
        printf("query number:%d, enough~\n", i);
        for(int k = 0; k < queries.size(); k++) {
            stringstream s_stream;
            unordered_set<int> sample = queries[k];
            for(const auto & element: sample) {
                s_stream << element << " ";
            }
            string s = s_stream.str();
            query_f<<s.substr(0, s.size() - 1)<<"\n";
        }
        query_f.close();
    }

    igraph_vector_destroy(&degrees);
}

void generate_queries_vary_distance(igraph_t* graph, string & query_dir) {
    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    vector<int> degree_arrays;
    vector<int> degree_starts;
    vector<int> degree_ends;
    get_sorted_degrees(graph, degrees, degree_arrays, degree_starts, degree_ends);

    dynamic_bitset<> candidates((unsigned long)igraph_vcount(graph));
    int start_index = (igraph_vector_size(&degrees) - degree_starts[1]) * 0.2 + degree_starts[1];
    int degree_start = degree_starts[VECTOR(degrees)[degree_arrays[start_index]]];

    int end_index = (igraph_vector_size(&degrees) - degree_starts[1]) * 0.4 + degree_starts[1];
    int degree_end = degree_ends[VECTOR(degrees)[degree_arrays[end_index]]];

    int size = igraph_vcount(graph);

    for(int i = degree_start; i <= degree_end; i++)
        candidates.set((unsigned long)(degree_arrays[i]));

    printf("candidates size:%d\n", (int)candidates.count());

    for(int i = 1; i <= 5; i++) {

        stringstream dir_stram;
        dir_stram<<query_dir<<"/vary_dist.query."<<i;

        fstream query_f(dir_stram.str(), ios::out);

        vector<unordered_set<int>> queries;
        while(queries.size() < 100) {
            unordered_set<int> query = vary_q_number(graph, i, 3, candidates);
            if(query.size() > 0 && !contains(queries, query)) {
                queries.push_back(query);
                printf("query size:%d\n", (int)queries.size());
            }
        }
        printf("query dist:%d, enough~\n", i);
        for(int k= 0; k < queries.size(); k++) {
            stringstream s_stream;
            unordered_set<int> sample = queries[k];
            for(const auto & element: sample) {
                s_stream<< element<< " ";
            }
            string s = s_stream.str();
            query_f<< s.substr(0, s.size() - 1)<<"\n";
        }
        query_f.close();
    }

    igraph_vector_destroy(&degrees);
}


void generate_default_queries(igraph_t* graph, string & query_dir) {

    fstream query_f((query_dir + "/default.query"), ios::out);

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    vector<int> degree_arrays;
    vector<int> degree_starts;
    vector<int> degree_ends;

    get_sorted_degrees(graph, degrees, degree_arrays, degree_starts, degree_ends);

    dynamic_bitset<> candidates((unsigned long)igraph_vcount(graph));
    int index = (igraph_vector_size(&degrees) - degree_starts[1]) * 0.2 + degree_starts[1];

    int degree_start = degree_starts[VECTOR(degrees)[degree_arrays[index]]];
    int degree_end = degree_ends[VECTOR(degrees)[degree_arrays[index]]];

    for(int i = degree_start; i <= degree_end; i++) {
        candidates.set((unsigned long) (degree_arrays[i]));
    }

    printf("candidates size:%d\n", (int)candidates.count());

    vector<unordered_set<int>> queries;
    while(queries.size() < 100) {
        unordered_set<int> query = vary_q_number(graph, 2, 3, candidates);
        if(query.size() > 0 && !contains(queries, query)) {
            queries.push_back(query);
            printf("query size:%d\n", (int)queries.size());
        }
    }
    printf("query, enough~\n");
    for(int k = 0; k < queries.size(); k++) {
        stringstream s_stream;
        unordered_set<int> sample = queries[k];
        for(const auto & element: sample) {
            s_stream << element << " ";
        }
        string s = s_stream.str();
        query_f<<s.substr(0, s.size() - 1)<<"\n";
    }

    igraph_vector_destroy(&degrees);

    query_f.close();
}



void generate_queries(string & file_name, string & query_directory) {

    srand((unsigned)time(NULL));

    igraph_t graph;

    get_simple_graph(file_name.data(), &graph);

    printf("generate default queires\n");
    generate_default_queries(&graph, query_directory);

    printf("generate varying degree queries\n");
    generate_queries_vary_degree(&graph, query_directory);

    printf("generate varying qnumber queries\n");
    generate_queries_vary_qnumber(&graph, query_directory);

    printf("generate varying qdist queries\n");
    generate_queries_vary_distance(&graph, query_directory);

    igraph_destroy(&graph);
}
