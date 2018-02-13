//
// Created by wz on 18-2-7.
//

#include "LocKDR.h"

igraph_t* loc_graph;
double epsilon = 0.03;

void init_loc(igraph_t* graph, int steiner_type) {
    loc_graph = graph;
    init_steiner(graph, steiner_type);
}

void set_c(double c) {
    set_steiner_c(c);
}

void set_alpha(double alpha) {
    set_steiner_alpha(alpha);
}


vector<vector<int>> compute_query_distances(igraph_t* graph,
                                            int & max_distance,
                                            unordered_map<int,int> & inverted_queires,
                                            unordered_map<int, int> & lookup_table) {
    max_distance = 0;

    vector<vector<int>> query_distances(igraph_vcount(graph));// vertex-> <query distance>, max_query_distance
    for(int i = 0; i < igraph_vcount(graph); i++) {
        query_distances[i].resize(inverted_queires.size() + 1);
        query_distances[i][inverted_queires.size()] = 0;
    }

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    for(const auto & query_pair: inverted_queires) {
        auto element = lookup_table[query_pair.first];
        query_distances[element][query_pair.second] = 0;

        queue<int> bfs_queue;
        bfs_queue.push(element);

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            auto source = bfs_queue.front();
            bfs_queue.pop();

            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            int distance = query_distances[source][query_pair.second] + 1;

            if(query_distances[source][inverted_queires.size()] < distance - 1)
                query_distances[source][inverted_queires.size()] = distance - 1;

            if(max_distance < distance -1)
                max_distance = distance - 1;

            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                    query_distances[neighbor][query_pair.second] = distance;
                }
            }
        }
    }

    igraph_vector_destroy(&neighbors);
    return query_distances;
}


void expand_steiner_vertices(unordered_set<int> & vertices,
                             vector<vector<int>> & query_distances,
                             unordered_map<int, int> & lookup_table,
                             unordered_map<int, int> & inverted_queries,
                             int upper_bound) {
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    // vertex-> <query_distance>, max_query_distance from vertices
    unordered_map<int, vector<int>> frontier_stats;
    for(const auto & vertex: vertices) {
        igraph_neighbors(loc_graph, &neighbors, vertex, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)VECTOR(neighbors)[i];
            if(vertices.find(neighbor) == vertices.end()) {
                if(frontier_stats.find(neighbor) == frontier_stats.end()) {
                    frontier_stats[neighbor].resize(inverted_queries.size() + 1);
                    int vertex_index = lookup_table[vertex];
                    for(int l = 0; l < inverted_queries.size(); l++)
                        frontier_stats[neighbor][l] = query_distances[vertex_index][l] + 1;
                    frontier_stats[neighbor][inverted_queries.size()]
                            = query_distances[vertex_index][inverted_queries.size()];
                } else {
                    int vertex_index = lookup_table[vertex];
                    for(int l = 0; l < inverted_queries.size(); l++) {
                        int query_distance = query_distances[vertex_index][l] + 1;
                        if(frontier_stats[neighbor][l] > query_distance)
                            frontier_stats[neighbor][l] = query_distance;
                    }
                    if(frontier_stats[neighbor][inverted_queries.size()]
                       < query_distances[vertex_index][inverted_queries.size()])
                        frontier_stats[neighbor][inverted_queries.size()]
                                = query_distances[vertex_index][inverted_queries.size()];
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);

    priority_queue<relNode, vector<relNode>, relNode_can_cmp> candidates;
    for(const auto & frontier_pair: frontier_stats) {
        int temp_max_distance = -1;
        for(int l = 0; l < inverted_queries.size(); l++) {
            if(temp_max_distance < frontier_pair.second[l])
                temp_max_distance = frontier_pair.second[l];
        }
        int super_max_distance = frontier_pair.second[inverted_queries.size()];
        if(temp_max_distance <= super_max_distance) {
            double relevance = -(temp_max_distance + 0.0) * temp_max_distance / super_max_distance;
            candidates.push(relNode(frontier_pair.first, relevance, 0));
        }
    }

    while(!candidates.empty() && vertices.size() < upper_bound) {
        vertices.insert(candidates.top().vertex);
        candidates.pop();
    }

}

unordered_set<int> get_removed_vertices(igraph_t* graph,
                                        unordered_set<int>* query_nodes,
                                        unordered_set<int>* removed_vertices,
                                        int & max_distance) {
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    int size = (int)(igraph_vcount(graph) - removed_vertices->size());
    double init_relevance = 1.0 / size;
    unordered_map<int, vector<double>> vertices_info;
    // preserved_vertex-> current_query_distance, max_query_distance,
    // older_relevance, newer_relevance, neighbor_size
    for(int i = 0; i < igraph_vcount(graph); i++) {
        if(removed_vertices->find(i) == removed_vertices->end()) {
            igraph_neighbors(graph, &neighbors, i, IGRAPH_ALL);
            int degree = 0;
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];
                if(removed_vertices->find(neighbor) == removed_vertices->end())
                    degree++;
            }
            vertices_info[i] = {0, 0, init_relevance, 0.0, (double)degree};
        }
    }

    int older_relevance_order = 2;
    int newer_relevance_order = 3;
    int swap_relevance_order = 2;
    for(const auto & element: *query_nodes) {
        queue<int> bfs_queue;
        bfs_queue.push(element);
        vertices_info[element][0] = 0;

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            bfs_queue.pop();
            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            double current_distance = vertices_info[source][0] + 1;
            vertices_info[source][newer_relevance_order] = 0.0;
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(removed_vertices->find(neighbor) == removed_vertices->end()) {
                    vertices_info[source][newer_relevance_order] +=
                            vertices_info[neighbor][older_relevance_order] / vertices_info[neighbor][4];
                    if(visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        if(current_distance > vertices_info[neighbor][1])
                            vertices_info[neighbor][1] = current_distance;
                        vertices_info[neighbor][0] = current_distance;
                        bfs_queue.push(neighbor);
                    }
                }
            }
        }

        swap_relevance_order = older_relevance_order;
        older_relevance_order = newer_relevance_order;
        newer_relevance_order = swap_relevance_order;
    }
    igraph_vector_destroy(&neighbors);

    priority_queue<relNode, vector<relNode>, relNode_bi_cmp> candidates;
    max_distance = -1;
    for(const auto & element: vertices_info) {
        if(max_distance <= element.second[1]) {
            max_distance = (int)element.second[1];
            candidates.push(relNode(element.first,
                                    element.second[older_relevance_order],
                                    (int)element.second[1]));
        }
    }

    size = (int)(epsilon * size);
    if(size <= 0)
        size = 1;

    unordered_set<int> vertices_to_del;
    while(vertices_to_del.size() < size && !candidates.empty()) {
        auto element = candidates.top();
        vertices_to_del.insert(element.vertex);
        candidates.pop();
    }

    return vertices_to_del;
}



unordered_map<int, vector<double>> compute_phps(igraph_t* graph,
                                                unordered_set<int> * query_nodes,
                                                unordered_set<int> * removed_vertices) {
    double php_epsilon = 0.01;
    double php_decay = 0.9;
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    unordered_map<int, vector<double>> vertices_info;

    int start_vertex = *(query_nodes->begin());

    unordered_set<int> visited;
    visited.insert(start_vertex);

    queue<int> bfs_queue;
    bfs_queue.push(start_vertex);

    while(!bfs_queue.empty()) {
        auto source = bfs_queue.front();
        bfs_queue.pop();

        double degree = 0;
        igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
        for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
            auto neighbor = (int)VECTOR(neighbors)[l];

            if(removed_vertices->find(neighbor)
               == removed_vertices->end()) {
                degree++;
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                }
            }
        }

        if(query_nodes->find(source) != query_nodes->end())
            vertices_info[source] = {1.0, 0.0, degree, 0, 0};//older_order, newer_order, degree, current_d, max_d
        else
            vertices_info[source] = {0.0, 0.0, degree, 0, 0};
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
        bfs_queue.push(start_vertex);
        visited.insert(start_vertex);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            bfs_queue.pop();

            vertices_info[source][newer_order] = 0.0;
            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];

                if(removed_vertices->find(neighbor)
                   == removed_vertices->end()) {
                    vertices_info[source][newer_order]
                            += vertices_info[neighbor][older_order];
                    if(visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        bfs_queue.push(neighbor);
                    }
                }
            }

            if(query_nodes->find(source) != query_nodes->end())
                vertices_info[source][newer_order] = 1.0;
            else {
                vertices_info[source][newer_order] =
                        vertices_info[source][newer_order] * php_decay / vertices_info[source][2];
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

    for(auto & pair: vertices_info)
        pair.second[0] = pair.second[older_order];

    igraph_vector_destroy(&neighbors);
    return vertices_info;
};


unordered_set<int> get_removed_vertices_php(igraph_t* graph, double & php_metric,
                                            unordered_set<int>* query_nodes,
                                            unordered_set<int>* removed_vertices,
                                            int & max_distance) {

    //older_order, newer_order, degree, current_d, max_d
    auto phps = compute_phps(graph, query_nodes, removed_vertices);

    php_metric = 0.0;
    for(const auto & pair:phps)
        php_metric += pair.second[0] * pair.second[2];
    php_metric = php_metric/phps.size()/(phps.size() - 1);


    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    for(const auto & element: *query_nodes) {
        queue<int> bfs_queue;
        bfs_queue.push(element);
        phps[element][3] = 0;

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            bfs_queue.pop();
            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            double current_distance = phps[source][3] + 1;
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(removed_vertices->find(neighbor) == removed_vertices->end()) {
                    if(visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        if(current_distance > phps[neighbor][4])
                            phps[neighbor][4] = current_distance;
                        phps[neighbor][3] = current_distance;
                        bfs_queue.push(neighbor);
                    }
                }
            }
        }
    }

    igraph_vector_destroy(&neighbors);


    priority_queue<relNode, vector<relNode>, relNode_bi_cmp> candidates;
    max_distance = -1;
    for(const auto & element: phps) {
        if(max_distance <= element.second[4]) {
            max_distance = (int)element.second[4];
            candidates.push(relNode(element.first,
                                    element.second[0] * element.second[2],
                                    (int)element.second[4]));
        }
    }

    int size = (int)(epsilon * phps.size());
    if(size <= 0)
        size = 1;

    unordered_set<int> vertices_to_del;
    while(vertices_to_del.size() < size && !candidates.empty()) {
        auto element = candidates.top();
        vertices_to_del.insert(element.vertex);
        candidates.pop();
    }

    return vertices_to_del;
}


void enlarge_target(double & kd_metric, double & php_metric,
                    unordered_set<int> & vertices,
                    unordered_set<int> & query_nodes,
                    igraph_vector_int_t* query_vertices) {
    igraph_t expanded_graph;
    vector<int> expanded_names;
    auto lookup_table = induce_subgraph_from_ori(loc_graph, &vertices,
                                                 &expanded_graph, &expanded_names);

    unordered_set<int> query_nodes_in_exp;
    for(const auto & element: query_nodes)
        query_nodes_in_exp.insert(lookup_table[element]);

    int maximal_core;
    igraph_t maximal_core_graph;
    vector<int> maximal_core_names;
    unordered_set<int> query_nodes_in_maximal_core;
    compute_maximal_core(&query_nodes_in_exp, &maximal_core,
                         &expanded_graph, &expanded_names,
                         &maximal_core_graph, &maximal_core_names,
                         query_nodes_in_maximal_core);
    igraph_destroy(&expanded_graph);

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(&maximal_core_graph, &degrees,
                  igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);

    int last_distance = INT32_MAX;
    double last_php_metric = -1.0;
    unordered_set<int> removed_vertices;
    unordered_set<int> last_removed_vertices;

    while(true) {
        int query_distance;
        double temp_php_metric;
        auto pending_vertices = get_removed_vertices_php(&maximal_core_graph, temp_php_metric,
                                                         &query_nodes_in_maximal_core,
                                                         &removed_vertices, query_distance);

        if(last_distance > query_distance ||
           (last_distance == query_distance && last_php_metric <= temp_php_metric)) {
            last_distance = query_distance;
            last_removed_vertices = removed_vertices;
            last_php_metric = temp_php_metric;
        }

        pending_vertices = maintain_k_core(&maximal_core_graph,
                                           &pending_vertices, &removed_vertices,
                                           &degrees, maximal_core);

        auto connected = check_query_connectivity(&maximal_core_graph,
                                                  *(query_nodes_in_maximal_core.begin()),
                                                  &query_nodes_in_maximal_core,
                                                  &removed_vertices, &pending_vertices);

        if(!connected)
            break;

        for(const auto & element: pending_vertices)
            removed_vertices.insert(element);
    }

    vertices.clear();
    int start_vertex = *(query_nodes_in_maximal_core.begin());
    auto query_component = get_query_component(&maximal_core_graph, start_vertex, &last_removed_vertices);


    for(const auto & element: query_component)
        vertices.insert(maximal_core_names[element]);

    igraph_vector_destroy(&degrees);
    igraph_destroy(&maximal_core_graph);

    kd_metric = (maximal_core + 0.0) / last_distance;
    php_metric = last_php_metric;
}




void expand_surrounding_vertices(int max_distance,
                                 double current_kd,
                                 int upper_bound,
                                 unordered_set<int> & vertices,
                                 vector<vector<int>> & query_distances,
                                 unordered_map<int, int> & inverted_queries,
                                 unordered_map<int, int> & lookup_table) {
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    unordered_map<int, vector<int>> frontier_stats;
    // frontier-> <query_distance> + k_degree
    for(const auto & vertex: vertices) {
        igraph_neighbors(loc_graph, &neighbors, vertex, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)VECTOR(neighbors)[i];
            if(vertices.find(neighbor) == vertices.end()) {
                if(frontier_stats.find(neighbor) == frontier_stats.end()) {
                    frontier_stats[neighbor].resize(inverted_queries.size() + 1);
                    int vertex_index = lookup_table[vertex];
                    for(int l = 0; l < inverted_queries.size(); l++) {
                        int query_distance = query_distances[vertex_index][l] + 1;
                        frontier_stats[neighbor][l] = query_distance;
                    }
                    frontier_stats[neighbor][inverted_queries.size()] = 1;
                } else {
                    int vertex_index = lookup_table[vertex];
                    for(int l = 0; l < inverted_queries.size(); l++) {
                        int query_distance = query_distances[vertex_index][l] + 1;
                        if(frontier_stats[neighbor][l] > query_distance)
                            frontier_stats[neighbor][l] = query_distance;
                    }
                    frontier_stats[neighbor][inverted_queries.size()]++;
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);

    priority_queue<relNode, vector<relNode>, relNode_can_cmp> candidates;
    for(const auto & frontier_pair: frontier_stats) {
        int temp_max_distance = -1;
        for(int l = 0; l < inverted_queries.size(); l++)
            if(temp_max_distance < frontier_pair.second[l])
                temp_max_distance = frontier_pair.second[l];
        if(temp_max_distance < max_distance)
            temp_max_distance = max_distance;
        double temp_k = frontier_pair.second[inverted_queries.size()] + 0.0;
        if(temp_k / temp_max_distance >= current_kd)
            candidates.push(relNode(frontier_pair.first, temp_k / temp_max_distance, 0));
    }

    while(!candidates.empty() && vertices.size() < upper_bound) {
        vertices.insert(candidates.top().vertex);
        candidates.pop();
    }
}


unordered_set<int> generate_loc_community(int upper_bound,
                                          unordered_set<int> & vertices,
                                          igraph_vector_int_t* query_vertices) {
    unordered_set<int> query_nodes;
    unordered_map<int, int> inverted_queries;// query_vertex-> index in array
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        query_nodes.insert(VECTOR(*query_vertices)[i]);
        inverted_queries[VECTOR(*query_vertices)[i]] = i;
    }

    // expand steiner tree
    igraph_t steiner_graph;
    vector<int> steiner_names;
    auto lookup_table = induce_subgraph_from_ori(loc_graph,
                                                 &vertices, &steiner_graph, &steiner_names);

    int max_distance;
    vector<vector<int>> query_distances =
            compute_query_distances(
                    &steiner_graph, max_distance, inverted_queries, lookup_table);

    expand_steiner_vertices(vertices, query_distances, lookup_table, inverted_queries, upper_bound);

    lookup_table.clear();
    igraph_destroy(&steiner_graph);

    // enlarge kd metric
    double last_metric;
    double last_php_metric;
    if(vertices.size() >= upper_bound) {
        enlarge_target(last_metric, last_php_metric, vertices, query_nodes, query_vertices);
        return vertices;
    }

    enlarge_target(last_metric, last_php_metric, vertices, query_nodes, query_vertices);
    auto last_vertices = vertices;
    double metric = last_metric;
    double php_metric = last_php_metric;

    // iteratively enlarging metric
    while(true) {
        igraph_t temp_graph;
        vector<int> temp_names;
        lookup_table = induce_subgraph_from_ori(loc_graph, &vertices, &temp_graph, &temp_names);
        query_distances = compute_query_distances(&temp_graph, max_distance, inverted_queries, lookup_table);
        expand_surrounding_vertices(max_distance, metric, upper_bound,
                                    vertices, query_distances,
                                    inverted_queries, lookup_table);
        igraph_destroy(&temp_graph);

        if(vertices.size() >= upper_bound) {
            enlarge_target(metric, php_metric, vertices, query_nodes, query_vertices);
            if(last_metric > metric
               || last_metric == metric && last_php_metric >= php_metric)
                return last_vertices;
            else
                return vertices;
        }

        enlarge_target(metric, php_metric, vertices, query_nodes, query_vertices);
        if(vertices.size() > upper_bound || last_metric > metric)
            break;

        if(last_metric == metric && last_php_metric >= php_metric)
            break;

        last_vertices = vertices;
        last_metric = metric;
        last_php_metric = php_metric;
    }
    return last_vertices;
}

vector<unordered_set<int>> generate_seeds(int vertex) {

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    igraph_neighbors(loc_graph, &neighbors, vertex, IGRAPH_ALL);

    unordered_set<int> vertices;
    for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
        auto neighbor = (int)VECTOR(neighbors)[i];
        vertices.insert(neighbor);
    }
    vertices.insert(vertex);

    igraph_t ego_graph;
    vector<int> ego_names;
    auto lookup_table = induce_subgraph_from_ori(
            loc_graph, &vertices, &ego_graph, &ego_names);
    vertices.erase(vertex);

    vector<unordered_set<int>> partitions;
    while(!vertices.empty()) {

        unordered_set<int> visited;
        int start_vertex = lookup_table[*(vertices.begin())];
        visited.insert(start_vertex);
        visited.insert(lookup_table[vertex]);
        queue<int> bfs_queue;
        bfs_queue.push(start_vertex);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            bfs_queue.pop();

            igraph_neighbors(&ego_graph, &neighbors, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                }
            }
        }

        unordered_set<int> partition;
        for(const auto & element: visited) {
            partition.insert(ego_names[element]);
            vertices.erase(ego_names[element]);
        }

        partitions.push_back(partition);
    }

    igraph_destroy(&ego_graph);
    igraph_vector_destroy(&neighbors);
    return partitions;
}

vector<unordered_set<int>> get_loc_kdr_single(
        igraph_vector_int_t* query_vertices, int upper_bound) {

    auto partitions = generate_seeds(VECTOR(*query_vertices)[0]);

    vector<unordered_set<int>> communities;
    for(int l = 0; l < partitions.size(); l++) {
        communities.push_back(generate_loc_community(upper_bound, partitions[l], query_vertices));
    }
    return communities;
}

unordered_set<int> get_loc_kdr(igraph_vector_int_t* query_vertices, int upper_bound) {

    // compute a steiner tree containing Q
    auto steiner_edges = get_steiner_tree(query_vertices);

    unordered_set<int> steiner_vertices;
    for(const auto & steiner_edge: steiner_edges) {
        int source, target;
        igraph_edge(loc_graph, steiner_edge, &source, &target);
        steiner_vertices.insert(source);
        steiner_vertices.insert(target);
    }

    return generate_loc_community(upper_bound, steiner_vertices, query_vertices);
}