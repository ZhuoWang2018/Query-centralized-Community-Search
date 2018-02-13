//
// Created by wz on 18-2-7.
//

#include "k_core_decomp.h"

void del_degree_vertex(int* degree_starts, igraph_vector_t* degrees, int vertex_idx) {
    auto degree = (int)(VECTOR(*degrees)[vertex_idx]);
    degree_starts[degree]++;
    VECTOR(*degrees)[vertex_idx]--;
}

void update_degree_vertex(int* degree_starts, igraph_vector_t* degrees,
                          int* degree_positions, int* degree_arrays, int vertex_idx) {
    auto degree = (int)(VECTOR(*degrees)[vertex_idx]);
    auto position_1 = degree_starts[degree];
    auto position_2 = degree_positions[vertex_idx];
    if(position_1 != position_2) {
        auto temp = degree_arrays[position_1];
        degree_arrays[position_1] = degree_arrays[position_2];
        degree_arrays[position_2] = temp;
        degree_positions[vertex_idx] = position_1;
        degree_positions[degree_arrays[position_2]] = position_2;
    }
    degree_starts[degree] += 1;
    VECTOR(*degrees)[vertex_idx]--;
}

bool check_query_connectivity(igraph_t* g, int start_vertex, unordered_set<int>* query_nodes,
                              unordered_set<int>* removed_vertices, unordered_set<int>* pending_vertices) {
    if(pending_vertices->find(start_vertex) != pending_vertices->end())
        return false;

    queue<int> bfs_queue;
    unordered_set<int> visited;
    int query_traveled_counter = 1;
    visited.insert(start_vertex);
    bfs_queue.push(start_vertex);

    if(query_traveled_counter >= query_nodes->size())
        return true;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(!bfs_queue.empty()) {
        auto source = bfs_queue.front();
        bfs_queue.pop();
        igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            int neighbor = (int)VECTOR(neighbors)[i];
            if(visited.find(neighbor) == visited.end() &&
               removed_vertices->find(neighbor) == removed_vertices->end() &&
               pending_vertices->find(neighbor) == pending_vertices->end()) {
                visited.insert(neighbor);
                bfs_queue.push(neighbor);
                if(query_nodes->find(neighbor) != query_nodes->end())
                    query_traveled_counter++;
                if(query_traveled_counter >= query_nodes->size()) {
                    igraph_vector_destroy(&neighbors);
                    return true;
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return false;
}

unordered_set<int> get_query_component(igraph_t* g, int start_vertex,
                                       unordered_set<int>* removed_vertices) {
    queue<int> vital_queue;
    unordered_set<int> visited;
    visited.insert(start_vertex);
    vital_queue.push(start_vertex);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(!vital_queue.empty()) {
        auto source = vital_queue.front();
        vital_queue.pop();
        igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)(VECTOR(neighbors)[i]);
            if(removed_vertices->find(neighbor) == removed_vertices->end() &&
               visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                vital_queue.push(neighbor);
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return visited;
}


int get_induced_vertex(int vertex_name, unordered_map<int, int> * induced_map) {
    if(induced_map->find(vertex_name) == induced_map->end()) {
        auto dict_size = induced_map->size();
        induced_map->insert({vertex_name, dict_size});
        return (int)dict_size;
    }
    return induced_map->at(vertex_name);
}

// vertex_in_ori:vertex_in_induced
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

unordered_map<int, int> induce_subgraph(
        unordered_set<int>* vertices,
        igraph_t* graph, vector<int>* names,
        igraph_t* induced_graph, vector<int>* induced_names) {

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);

    unordered_set<int> visited_edges;

    unordered_map<int, int> result;

    for(const auto & element: *vertices) {
        igraph_incident(graph, &neighbor_eids, element, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
            if(visited_edges.find(neighbor_eid) == visited_edges.end()) {
                int from = IGRAPH_FROM(graph, neighbor_eid);
                int to = IGRAPH_TO(graph, neighbor_eid);
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
        induced_names->at(element.second) = names->at(element.first);

    return result;
};


bool compute_maximal_k_core_from_ori(unordered_set<int> & query_nodes,
                                     igraph_vector_int_t* query_vertices,
                                     int k_core, igraph_t* g,
                                     igraph_t* k_core_graph, vector<int>* k_core_names,
                                     unordered_set<int> & query_nodes_in_k_core) {
    int bfs_start_vertex = VECTOR(*query_vertices)[0];

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    auto degree_arrays = new int[igraph_vcount(g)]();

    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    auto degree_starts = new int[max_degree + 1]();
    auto degree_ends = new int[max_degree + 1]();
    auto vertex_positions = new int[igraph_vcount(g)]();

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(g); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        vertex_positions[i] = degree_ends[vertex_degree];
        degree_ends[vertex_degree]++;
    }
    delete degree_ends;

    unordered_set<int> removed_vertices;
    unordered_set<int> pending_vertices;

    auto k = (int)(VECTOR(degrees)[degree_arrays[0]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(g);
    bool k_exceed_error = false;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(k < k_core) {
        if(live_start >= len_degree_arrays) {
            k_exceed_error = true;
            break;
        }

        int vertex_idx = degree_arrays[live_start];

        if((int)(VECTOR(degrees)[vertex_idx]) <= k) {
            del_degree_vertex(degree_starts, &degrees, vertex_idx);
            igraph_neighbors(g, &neighbors, vertex_idx, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)(VECTOR(neighbors)[i]);
                if((int)(VECTOR(degrees)[neighbor]) > k) {
                    update_degree_vertex(degree_starts, &degrees, vertex_positions, degree_arrays, neighbor);
                    if((int)(VECTOR(degrees)[neighbor]) <= k &&
                       query_nodes.find(neighbor) != query_nodes.end()) {
                        k_exceed_error = true;
                        break;
                    }
                }
            }
            pending_vertices.insert(vertex_idx);
            if(k_exceed_error)
                break;
            live_start = degree_starts[k];
        } else {
            bool checked = check_query_connectivity(g, bfs_start_vertex,
                                                    &query_nodes, &removed_vertices, &pending_vertices);

            if(checked) {
                for(const auto & element: pending_vertices)
                    removed_vertices.insert(element);
                pending_vertices.clear();
            } else {
                k_exceed_error = true;
                break;
            }
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }

    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

    if(k_exceed_error)
        return false;

    auto vertices = get_query_component(g, bfs_start_vertex, &removed_vertices);
    auto lookup_tales = induce_subgraph_from_ori(g, &vertices, k_core_graph, k_core_names);
    for(const auto & element: query_nodes)
        query_nodes_in_k_core.insert(lookup_tales[element]);

    return true;
}


void k_core_decomposition(igraph_t* g) {

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(g, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    auto degree_arrays = new int[igraph_vcount(g)]();

    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    auto degree_starts = new int[max_degree + 1]();
    auto degree_ends = new int[max_degree + 1]();
    auto vertex_positions = new int[igraph_vcount(g)]();

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(g); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        vertex_positions[i] = degree_ends[vertex_degree];
        degree_ends[vertex_degree]++;
    }
    delete degree_ends;

    auto k = (int)(VECTOR(degrees)[degree_arrays[0]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(g);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(true) {
        if(live_start % 100000 == 0 && live_start > 0)
            printf("I am still computing cores:%d\n", live_start);

        if(live_start >= len_degree_arrays)
            break;
        int vertex_idx = degree_arrays[live_start];
        if((int)(VECTOR(degrees)[vertex_idx]) <= k) {
            del_degree_vertex(degree_starts, &degrees, vertex_idx);
            igraph_neighbors(g, &neighbors, vertex_idx, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if((int)(VECTOR(degrees)[neighbor]) > k)
                    update_degree_vertex(degree_starts, &degrees,
                                         vertex_positions, degree_arrays, neighbor);
            }
            live_start = degree_starts[k];
            SETVAN(g, "cnum", vertex_idx, k);
        } else
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
    }

    SETGAN(g, "max-cnum", k);
    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

}


void compose_maximal_core(igraph_t* g,
                          int & cur_cnum, int & next_cnum,
                          unordered_set<int> & query_nodes,
                          igraph_vector_int_t* query_vertices,
                          igraph_t* k_core_graph, vector<int>* k_core_names,
                          unordered_set<int> & query_nodes_in_k_core) {

    int maximal_k_core = INT32_MAX;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++)
        if (maximal_k_core > VAN(g, "cnum", VECTOR(*query_vertices)[i]))
            maximal_k_core = (int) (VAN(g, "cnum", VECTOR(*query_vertices)[i]));

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    int travel_counter;
    unordered_set<int> vertices;

    do {
        vertices.clear();

        travel_counter = 1;
        next_cnum = -1;

        int start_vertex = VECTOR(*query_vertices)[0];
        queue<int> bfs_queue;
        unordered_set<int> visited;
        bfs_queue.push(start_vertex);
        visited.insert(start_vertex);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front();
            vertices.insert(source);
            bfs_queue.pop();

            igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];
                if(visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    auto current_cnum = (int)VAN(g, "cnum", neighbor);
                    if(current_cnum >= maximal_k_core) {
                        bfs_queue.push(neighbor);
                        if(query_nodes.find(neighbor) != query_nodes.end())
                            travel_counter++;
                    } else if(current_cnum > next_cnum)
                        next_cnum = current_cnum;
                }
            }
        }
        cur_cnum = maximal_k_core;
        maximal_k_core = next_cnum;
    } while (travel_counter < query_nodes.size());

    igraph_vector_destroy(&neighbors);

    auto lookup_tables = induce_subgraph_from_ori(g, &vertices, k_core_graph, k_core_names);
    for(const auto & element: query_nodes)
        query_nodes_in_k_core.insert(lookup_tables[element]);
}

void compose_k_core(igraph_t* g, int cnum, int & next_cnum,
                    igraph_vector_int_t* query_vertices,
                    igraph_t* k_core_graph, vector<int>* k_core_names,
                    unordered_set<int> & query_nodes_in_k_core) {

    unordered_set<int> vertices;
    next_cnum = -1;

    int start_vertex = VECTOR(*query_vertices)[0];
    queue<int> bfs_queue;
    unordered_set<int> visited;
    bfs_queue.push(start_vertex);
    visited.insert(start_vertex);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(!bfs_queue.empty()) {
        int source = bfs_queue.front();
        vertices.insert(source);
        bfs_queue.pop();

        igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
        for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
            auto neighbor = (int)VECTOR(neighbors)[l];
            if(visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                auto current_cnum = (int)VAN(g, "cnum", neighbor);
                if(current_cnum >= cnum)
                    bfs_queue.push(neighbor);
                else if(current_cnum > next_cnum)
                    next_cnum = current_cnum;
            }
        }
    }

    igraph_vector_destroy(&neighbors);
    auto lookup_tables = induce_subgraph_from_ori(g, &vertices, k_core_graph, k_core_names);
    for(int l = 0; l < igraph_vector_int_size(query_vertices); l++)
        query_nodes_in_k_core.insert(lookup_tables[VECTOR(*query_vertices)[l]]);
};


void compute_maximal_core(
        unordered_set<int>* query_nodes_in_cur,
        int* maximal_core, igraph_t* graph, vector<int>* names,
        igraph_t* maximal_core_graph, vector<int>* maximal_core_names,
        unordered_set<int> & query_nodes_in_maximal_core) {

    int start_vertex = *(query_nodes_in_cur->begin());

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    auto degree_arrays = new int[igraph_vcount(graph)]();

    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    auto degree_starts = new int[max_degree + 1]();
    auto degree_ends = new int[max_degree + 1]();
    auto vertex_positions = new int[igraph_vcount(graph)]();

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(graph); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        vertex_positions[i] = degree_ends[vertex_degree];
        degree_ends[vertex_degree]++;
    }
    delete degree_ends;

    unordered_set<int> removed_vertices;
    unordered_set<int> pending_vertices;

    auto k = (int)(VECTOR(degrees)[degree_arrays[0]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(graph);
    bool turn_to_outer = false;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(true) {
        if(live_start >= len_degree_arrays)
            break;
        int vertex_idx = degree_arrays[live_start];
        if((int)(VECTOR(degrees)[vertex_idx]) <= k) {
            del_degree_vertex(degree_starts, &degrees, vertex_idx);
            igraph_neighbors(graph, &neighbors, vertex_idx, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)(VECTOR(neighbors)[i]);
                if((int)(VECTOR(degrees)[neighbor]) > k) {
                    update_degree_vertex(degree_starts, &degrees, vertex_positions, degree_arrays, neighbor);
                    if((int)(VECTOR(degrees)[neighbor]) <= k &&
                       query_nodes_in_cur->find(neighbor) != query_nodes_in_cur->end()) {
                        turn_to_outer = true;
                        break;
                    }
                }
            }
            pending_vertices.insert(vertex_idx);
            if(turn_to_outer)
                break;
            live_start = degree_starts[k];
        } else {
            bool checked = check_query_connectivity(graph, start_vertex,
                                                    query_nodes_in_cur,
                                                    &removed_vertices, &pending_vertices);

            if(checked) {
                for(const auto & element: pending_vertices)
                    removed_vertices.insert(element);
                pending_vertices.clear();
            } else
                break;
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }


    auto vital_vertices = get_query_component(graph, start_vertex, &removed_vertices);

    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

    auto lookup_table = induce_subgraph(&vital_vertices, graph, names,
                                        maximal_core_graph, maximal_core_names);

    for(const auto & element: *query_nodes_in_cur)
        query_nodes_in_maximal_core.insert(lookup_table[element]);


    *maximal_core = k;
}
