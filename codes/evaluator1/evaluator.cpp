//
// Created by wz on 18-1-16.
//
#include "evaluator.h"

igraph_t * graph;

void init(igraph_t* g) {
    graph = g;
}

void split(string & str, unordered_set<int> & vertices) {
    string::size_type pos;
    vertices.clear();
    string pattern = " ";
    str += pattern;
    unsigned long size = (int)str.size();

    for(unsigned long i = 0; i < size; i++) {
        pos = str.find(" ", i);
        if(pos < size) {
            string s = str.substr(i, pos - i);
            int number;
            sscanf(s.data(), "%d", &number);
            vertices.insert(number);
            i = pos + pattern.size() - 1;
        }
    }
}

int get_induced_vertex(int vertex_name, unordered_map<int, int> * induced_map) {
    if(induced_map->find(vertex_name) == induced_map->end()) {
        auto dict_size = induced_map->size();
        induced_map->insert({vertex_name, dict_size});
        return (int)dict_size;
    }
    return induced_map->at(vertex_name);
}


unordered_map<int, int> induce_subgraph_from_ori(igraph_t* g, unordered_set<int>* vertices,
                                                 igraph_t* induced_graph, vector<int>* induced_names) {
    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);

    unordered_set<int> visited_edges;

    unordered_map<int, int> result;

    for(const auto & element: *vertices) {
        igraph_incident(g, &neighbor_eids, element, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
            if(visited_edges.find(neighbor_eid) == visited_edges.end()) {
                int from = IGRAPH_FROM(g, neighbor_eid);
                int to = IGRAPH_TO(g, neighbor_eid);
                if(vertices->find(from) != vertices->end()
                   && vertices->find(to) != vertices->end()) {
                    auto from_vertex_idx = get_induced_vertex(from, &result);
                    auto to_vertex_idx = get_induced_vertex(to, &result);
                    igraph_vector_push_back(&edges, from_vertex_idx);
                    igraph_vector_push_back(&edges, to_vertex_idx);
                    visited_edges.insert(neighbor_eid);
                }
            }
        }
    }

    igraph_create(induced_graph, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&neighbor_eids);
    igraph_vector_destroy(&edges);

    induced_names->resize(vertices->size());
    for(const auto & element: result)
        induced_names->at((unsigned long)element.second) = element.first;

    return result;
}

int compute_max_distance(igraph_t* g, unordered_set<int>* nodes_in_g) {
    int max_distance = -1;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    for(const auto & element: *nodes_in_g) {
        queue<pair<int,int>> bfs_queue;
        bfs_queue.push(make_pair(element, 0));

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front().first;
            int distance = bfs_queue.front().second;
            bfs_queue.pop();

            if(distance > max_distance)
                max_distance = distance;

            igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(make_pair(neighbor, distance + 1));
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return max_distance;
}


void compute_metrics(string & s,
                     unordered_set<int> & vertices,
                     unordered_set<int> & query_nodes,
                     double & k, double & d, double & density,
                     double & metric_1, double & metric_2, double & metric_3) {
    split(s, vertices);

    igraph_t induced_graph;
    vector<int> induced_names;
    auto lookup_table = induce_subgraph_from_ori(
            graph, &vertices, &induced_graph, &induced_names);

    unordered_set<int> nodes_in_g;
    for(const auto & element: query_nodes)
        nodes_in_g.insert(lookup_table[element]);

    d = compute_max_distance(&induced_graph, &nodes_in_g);

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(&induced_graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    k = igraph_vector_min(&degrees);
    igraph_vector_destroy(&degrees);

    density = 2.0 * igraph_ecount(&induced_graph)
              / igraph_vcount(&induced_graph)
              / (igraph_vcount(&induced_graph) - 1);

    metric_1 = (k+0.0) / d;
    metric_2 = (k+0.0) * k / d;
    metric_3 = k + (k+0.0) / d;

    //printf("k:%lf, d:%lf, density:%lf, k/d:%lf, kk/d:%lf, k+k/d:%lf\n", k, d, density, metric_1, metric_2, metric_3);


    igraph_destroy(&induced_graph);
}




void handle_one_community(string & s,
                          unordered_set<int> & query_nodes,
                          unordered_set<int> & vertices,
                          double & sum_k, double & sum_d, double & sum_density,
                          double & sum_metric_1, double & sum_metric_2, double & sum_metric_3) {
    double k, d, density, metric_1, metric_2, metric_3;
    compute_metrics(s, vertices, query_nodes, k, d, density, metric_1, metric_2, metric_3);
    sum_k += k;
    sum_d += d;
    sum_density += density;
    sum_metric_1 += metric_1;
    sum_metric_2 += metric_2;
    sum_metric_3 += metric_3;
}

void handle_multi_communities(fstream & result_f, string & s, unordered_set<int> & query_nodes,
                              double & sum_k, double & sum_d, double & sum_density,
                              double & sum_metric_1, double & sum_metric_2, double & sum_metric_3) {

    double max_k = 0, max_d, max_density,
            max_metric_1 = 0, max_metric_2, max_metric_3;

    unordered_set<int> vertices;

    while(getline(result_f, s)) {
        if(s[0] != '#') {
            split(s, vertices);
            double k, d, density, metric_1, metric_2, metric_3;
            compute_metrics(s, vertices, query_nodes,
                            k, d, density, metric_1, metric_2, metric_3);
            if(max_metric_1 <= metric_1 && max_k <= k) {
                max_k = k;
                max_d = d;
                max_density = density;
                max_metric_1 = metric_1;
                max_metric_2 = metric_2;
                max_metric_3 = metric_3;
            }
        } else
            break;
    }

    sum_k += max_k;
    sum_d += max_d;
    sum_density += max_density;
    sum_metric_1 += max_metric_1;
    sum_metric_2 += max_metric_2;
    sum_metric_3 += max_metric_3;

}


void evaluate(string & query_file, string & result_file) {

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::in);

    double sum_k = 0.0;
    double sum_d = 0.0;
    double sum_density = 0.0;
    double sum_metric_1 = 0.0; // k/d
    double sum_metric_2 = 0.0; // k*k/d
    double sum_metric_3 = 0.0; // k + k/d

    string s;

    int counter = 0;

    while(getline(query_f, s)) {
        unordered_set<int> query_nodes;
        split(s, query_nodes);

        getline(result_f, s);
        if(s[0] == '#') {
            handle_multi_communities(result_f, s, query_nodes,
                sum_k, sum_d, sum_density,
                sum_metric_1, sum_metric_2, sum_metric_3);
        } else {
            unordered_set<int> vertices;
            split(s, vertices);
            handle_one_community(s, query_nodes, vertices,
                sum_k, sum_d, sum_density,
                sum_metric_1, sum_metric_2, sum_metric_3);
        }

        counter++;
    }

    query_f.close();
    result_f.close();

    double avg_k = sum_k / counter;
    double avg_dist = sum_d / counter;
    double avg_density = sum_density / counter;
    double avg_metric_1 = sum_metric_1 / counter;
    double avg_metric_2 = sum_metric_2 / counter;
    double avg_metric_3 = sum_metric_3 / counter;

    printf("avg-k:%lf, avg-dist:%lf, avg-k/d:%lf\n",
        avg_k, avg_dist, avg_metric_1);

}
