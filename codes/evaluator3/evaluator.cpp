//
// Created by wz on 18-1-30.
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

vector<double> compute_phps(igraph_t* graph,
                            unordered_set<int> * query_nodes,
                            igraph_vector_t & degrees) {
    double php_epsilon = 0.01;
    double php_decay = 0.9;
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

//    unordered_map<int, vector<double>> vertices_info;
    vector<vector<double>> vertices_info(igraph_vcount(graph));
    for(int i = 0; i < igraph_vcount(graph); i++) {
        if(query_nodes->find(i) != query_nodes->end())
            vertices_info[i] = {1.0, 0.0};
        else
            vertices_info[i] = {0.0, 0.0};
    }

    int older_order = 0;
    int newer_order = 1;
    int swap_order = 0;


    long passes = 10;
    if(passes < query_nodes->size())
        passes = query_nodes->size();

    int counter = 0;
    while(counter < passes) {
        double max_epsilon = -1.0;
        queue<int> bfs_queue;
        unordered_set<int> visited;
        bfs_queue.push(*(query_nodes->begin()));
        visited.insert(*(query_nodes->begin()));

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            bfs_queue.pop();

            vertices_info[source][newer_order] = 0.0;
            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];

                vertices_info[source][newer_order]
                        += vertices_info[neighbor][older_order];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                }
            }

            if(query_nodes->find(source) != query_nodes->end())
                vertices_info[source][newer_order] = 1.0;
            else {
                vertices_info[source][newer_order] =
                        vertices_info[source][newer_order] * php_decay / VECTOR(degrees)[source];
                if(max_epsilon < vertices_info[source][newer_order]
                                 - vertices_info[source][older_order])
                    max_epsilon = vertices_info[source][newer_order] - vertices_info[source][older_order];
            }
        }

        swap_order = newer_order;
        newer_order = older_order;
        older_order = swap_order;

        if(max_epsilon < php_epsilon)
            break;

        counter++;
    }


    vector<double> phps(vertices_info.size());
    for(int l = 0; l < vertices_info.size(); l++)
        phps[l] = vertices_info[l][older_order];

    igraph_vector_destroy(&neighbors);
    return phps;
};

double compute_weighted_density(
        unordered_set<int> & vertices,
        unordered_set<int> & query_nodes) {

    double weighted_density = 0.0;

    igraph_t induced_graph;
    vector<int> induced_names;
    auto lookup_table = induce_subgraph_from_ori(
            graph, &vertices, &induced_graph, &induced_names);

    unordered_set<int> nodes_in_g;
    for(const auto & element: query_nodes)
        nodes_in_g.insert(lookup_table[element]);

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);

    igraph_degree(&induced_graph, &degrees,
                  igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    auto phps = compute_phps(&induced_graph, &nodes_in_g, degrees);



    for(int l = 0; l < igraph_vector_size(&degrees); l++)
        weighted_density += (VECTOR(degrees)[l] * phps[l]);

    weighted_density = weighted_density / igraph_vcount(&induced_graph) / (igraph_vcount(&induced_graph) - 1);

    igraph_vector_destroy(&degrees);
    igraph_destroy(&induced_graph);
    return weighted_density;
}

double compute_max_weighted_density(fstream & result_f, string & s, unordered_set<int> & query_nodes) {
    double max_weighted_density = -1.0;

    unordered_set<int> vertices;

    while(getline(result_f, s)) {
        if(s[0] != '#') {
            split(s, vertices);
            auto weighted_density = compute_weighted_density(vertices, query_nodes);
            if(max_weighted_density < weighted_density)
                max_weighted_density = weighted_density;
        } else
            break;
    }

    return max_weighted_density;
}

void evaluate_weighted_density(string & query_file, string & result_file, bool single_included) {
    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::in);

    double weighted_density = 0.0;
    string s;

    int counter = 0;

    while(getline(query_f, s)) {
        unordered_set<int> query_nodes;
        split(s, query_nodes);

        getline(result_f, s);
        if(query_nodes.size() == 1) {
            if(s[0] == '#') {
                auto temp_weighted_density = compute_max_weighted_density(result_f, s, query_nodes);
                if(single_included) {
                    weighted_density += temp_weighted_density;
                    counter++;
                }
            } else {
                unordered_set<int> vertices;
                split(s, vertices);
                auto temp_weighted_density = compute_weighted_density(vertices, query_nodes);
                if(single_included) {
                    weighted_density += temp_weighted_density;
                    counter++;
                }
            }
        } else {
            unordered_set<int> vertices;
            split(s, vertices);
            weighted_density += compute_weighted_density(vertices, query_nodes);
            counter++;
        }
    }

    query_f.close();
    result_f.close();

    weighted_density /= counter;
    printf("weighted density:%lf\n", weighted_density);
}