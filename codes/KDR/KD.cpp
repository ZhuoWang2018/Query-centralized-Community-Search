//
// Created by wz on 18-2-7.
//
#include "KD.h"

igraph_t* kd_graph;

void init_kd(igraph_t* graph) {
    kd_graph = graph;
}

int get_query_distance_in_ori(igraph_vector_int_t* query_vertices) {

    int max_distance = -1;
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    for(int l = 0; l < igraph_vector_int_size(query_vertices) - 1; l++) {
        unordered_set<int> query_nodes;
        for(int m = l + 1; m < igraph_vector_int_size(query_vertices); m++)
            query_nodes.insert(VECTOR(*query_vertices)[m]);

        int start_vertex = VECTOR(*query_vertices)[l];
        unordered_set<int> visited;
        visited.insert(start_vertex);
        queue<pair<int,int>> bfs_queue;
        bfs_queue.push(make_pair(start_vertex, 0));

        int travel_counter = 0;
        while(!bfs_queue.empty()) {
            int source = bfs_queue.front().first;
            int distance = bfs_queue.front().second;
            bfs_queue.pop();

            int new_distance = distance + 1;
            igraph_neighbors(kd_graph, &neighbors, source, IGRAPH_ALL);
            for(int m = 0; m < igraph_vector_size(&neighbors); m++) {
                auto neighbor = (int)VECTOR(neighbors)[m];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(make_pair(neighbor, new_distance));
                    if(query_nodes.find(neighbor) != query_nodes.end()) {
                        travel_counter++;
                        if(new_distance > max_distance)
                            max_distance = new_distance;
                        if(travel_counter >= query_nodes.size())
                            break;

                    }
                }
            }

            if(travel_counter >= query_nodes.size())
                break;
        }
    }
    igraph_vector_destroy(&neighbors);
    return max_distance;
}


unordered_set<int> get_R_k(igraph_t* k_core_graph, vector<int>* k_core_names, int k_core, int & next_k_core,
                           unordered_set<int>* query_nodes_in_k_core, double & kd_metric, bool bulk) {
    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(k_core_graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    next_k_core = (int)igraph_vector_min(&degrees);

    int last_distance = INT32_MAX;
    unordered_set<int> removed_vertices;
    unordered_set<int> last_removed_vertices;

    while(true) {
        int query_distance;
        auto pending_vertices =
                get_max_query_dist_vertices(k_core_graph,
                                            query_nodes_in_k_core, &removed_vertices,
                                            query_distance, bulk);

        if (last_distance >= query_distance) {
            last_distance = query_distance;
            last_removed_vertices = removed_vertices;
        }


        pending_vertices = maintain_k_core(k_core_graph,
                                           &pending_vertices, &removed_vertices,
                                           &degrees, k_core);
        auto connected = check_query_connectivity(k_core_graph,
                                                  *(query_nodes_in_k_core->begin()),
                                                  query_nodes_in_k_core,
                                                  &removed_vertices, &pending_vertices);

        if (!connected)
            break;

        for(const auto & element: pending_vertices)
            removed_vertices.insert(element);

    }

    unordered_set<int> vertices;
    int start_vertex = *(query_nodes_in_k_core->begin());
    auto query_component = get_query_component(k_core_graph, start_vertex, &last_removed_vertices);

    int minimum_degree = INT32_MAX;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    for(const auto & element: query_component) {
        vertices.insert((k_core_names->at(element)));

        auto degree = 0;
        igraph_neighbors(k_core_graph, &neighbors, element, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)VECTOR(neighbors)[i];
            if(query_component.find(neighbor) != query_component.end())
                degree++;
        }

        if(degree < minimum_degree)
            minimum_degree = degree;

    }

    igraph_vector_destroy(&neighbors);
    igraph_vector_destroy(&degrees);

    kd_metric = (minimum_degree + 0.0) / last_distance;
    return vertices;

}



unordered_set<int> get_kd(igraph_vector_int_t* query_vertices) {
    unordered_set<int> query_nodes;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++)
        query_nodes.insert(VECTOR(*query_vertices)[i]);

    int cur_cnum;
    int next_cnum;
    igraph_t k_max_core_graph;
    vector<int> k_max_core_names;
    unordered_set<int> query_nodes_in_k_core;
    compose_maximal_core(kd_graph,
                         cur_cnum, next_cnum,
                         query_nodes, query_vertices,
                         &k_max_core_graph, &k_max_core_names,
                         query_nodes_in_k_core);

    double last_metric;
    int useless;
    unordered_set<int> last_vertices = get_R_k(&k_max_core_graph, &k_max_core_names,
                                               cur_cnum, useless,
                                               &query_nodes_in_k_core, last_metric, true);

    k_max_core_names.clear();
    igraph_destroy(&k_max_core_graph);

    int D_min;
    if(igraph_vector_int_size(query_vertices) == 1)
        D_min = 1;
    else
        D_min = get_query_distance_in_ori(query_vertices);

    cur_cnum = next_cnum;
    while(cur_cnum >= 1) {
        if((cur_cnum + 0.0) / D_min <= last_metric)
            break;

        double metric;
        igraph_t k_core_graph;
        vector<int> k_core_names;
        query_nodes_in_k_core.clear();
        compose_k_core(kd_graph,
                       cur_cnum, next_cnum,
                       query_vertices,
                       &k_core_graph, &k_core_names,
                       query_nodes_in_k_core);
        unordered_set<int> vertices = get_R_k(&k_core_graph, &k_core_names,
                                              cur_cnum, useless,
                                              &query_nodes_in_k_core, metric, true);
        igraph_destroy(&k_core_graph);

        if(last_metric < metric) {
            last_metric = metric;
            last_vertices = vertices;
        }
        cur_cnum = next_cnum;
    }

    return last_vertices;
}
