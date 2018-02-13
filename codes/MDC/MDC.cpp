//
// Created by wz on 18-1-17.
//
#include "MDC.h"

igraph_t copy_graph;
void init_mdc_ori(igraph_t * g_data) {
    igraph_copy(&copy_graph, g_data);
}

void destroy_mdc_ori() {
    igraph_destroy(&copy_graph);
}

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


bool bfs_connected_queries(igraph_t* graph,
                           int vertex_start_idx,
                           unordered_set<int>* vertex_idx_set,
                           boost::dynamic_bitset<>* removed_vertices,
                           boost::dynamic_bitset<>* pending_vertices) {
    if(pending_vertices->test((size_t)vertex_start_idx))
        return false;

    queue<int> bfs_queue;
    boost::dynamic_bitset<> visited(removed_vertices->size());
    int query_traveled_counter = 1;
    visited.set((size_t)vertex_start_idx);
    bfs_queue.push(vertex_start_idx);
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(!bfs_queue.empty()) {
        auto source = bfs_queue.front();
        bfs_queue.pop();
        igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            int neighbor = (int)VECTOR(neighbors)[i];
            if(!visited.test((size_t)neighbor)
               && !removed_vertices->test((size_t)neighbor)
               && !pending_vertices->test((size_t)neighbor)) {
                visited.set((size_t)neighbor);
                bfs_queue.push(neighbor);
                if(vertex_idx_set->find(neighbor) != vertex_idx_set->end())
                    query_traveled_counter++;
                if(query_traveled_counter >= vertex_idx_set->size()) {
                    igraph_vector_destroy(&neighbors);
                    return true;
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return false;
}

unordered_set<int> get_vital_component(igraph_t* graph,
                                       int vertex_start_idx,
                                       boost::dynamic_bitset<>* removed_vertices) {
    queue<int> vital_queue;
    unordered_set<int> visited;
    visited.insert(vertex_start_idx);
    vital_queue.push(vertex_start_idx);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(!vital_queue.empty()) {
        auto source = vital_queue.front();
        vital_queue.pop();
        igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)(VECTOR(neighbors)[i]);
            if(!removed_vertices->test(neighbor) &&
               visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                vital_queue.push(neighbor);
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return visited;
}

int get_induced_vertex_mdc(int vertex_name, unordered_map<int, int> * induced_map) {
    if(induced_map->find(vertex_name) == induced_map->end()) {
        auto dict_size = induced_map->size();
        induced_map->insert({vertex_name, dict_size});
        return (int)dict_size;
    }
    return induced_map->at(vertex_name);
}

unordered_set<int> induce_vertex_subgraph_mdc(igraph_t * graph, unordered_set<int>* filter_points, igraph_t* induced_graph,
                                              vector<int>* induced_names, unordered_set<int> & ori_query_set) {
    unordered_map<int, int> point_names_mdc;

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);

    unordered_set<int> visited_edges;

    for(const auto & element: *filter_points) {
        igraph_incident(graph, &neighbor_eids, element, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
            if(visited_edges.find(neighbor_eid) == visited_edges.end()) {
                int from = IGRAPH_FROM(graph, neighbor_eid);
                int to = IGRAPH_TO(graph, neighbor_eid);
                if(filter_points->find(from) != filter_points->end()
                   && filter_points->find(to) != filter_points->end()) {
                    auto from_vertex_idx = get_induced_vertex_mdc(from, &point_names_mdc);
                    auto to_vertex_idx = get_induced_vertex_mdc(to, &point_names_mdc);
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

    unordered_set<int> new_query_set;

    induced_names->resize(filter_points->size());
//    printf("induced names size:%ld\n", induced_names->size());
    for(const auto & element: point_names_mdc) {
        if(ori_query_set.find(element.first) != ori_query_set.end())
            new_query_set.insert(element.second);
        induced_names->at((unsigned long) element.second) = element.first;
    }
    return new_query_set;
}

// 0: have induced graph, induced names
// 1: use copied graph
int initial_decomp_max_core(igraph_vector_int_t * query_vertices, igraph_t* ori_graph,
                            igraph_t* induced_graph, vector<int>* induced_names, unordered_set<int> & vital_tool) {
    int start_idx = VECTOR(*query_vertices)[0];
    unordered_set<int> query_vertex_set;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        auto query_vertex = VECTOR(*query_vertices)[i];
        query_vertex_set.insert(query_vertex);
    }

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);
    igraph_degree(ori_graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    auto degree_arrays = new int[igraph_vcount(ori_graph)]();

    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    auto degree_starts = new int[max_degree + 1]();
    auto degree_ends = new int[max_degree + 1]();
    auto vertex_positions = new int[igraph_vcount(ori_graph)]();

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(ori_graph); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        vertex_positions[i] = degree_ends[vertex_degree];
        degree_ends[vertex_degree]++;
    }
    delete degree_ends;

    boost::dynamic_bitset<> pending_vertices((size_t)igraph_vcount(ori_graph));
    boost::dynamic_bitset<> removed_vertices((size_t)igraph_vcount(ori_graph));

    auto k = (int)(VECTOR(degrees)[degree_arrays[0]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(ori_graph);
    bool turn_to_outer = false;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    while(true) {
        if(live_start >= len_degree_arrays)
            break;
        int vertex_idx = degree_arrays[live_start];
        if((int)(VECTOR(degrees)[vertex_idx]) <= k) {
            del_degree_vertex(degree_starts, &degrees, vertex_idx);
            igraph_neighbors(ori_graph, &neighbors, vertex_idx, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)(VECTOR(neighbors)[i]);
                if((int)(VECTOR(degrees)[neighbor]) > k) {
                    update_degree_vertex(degree_starts, &degrees, vertex_positions, degree_arrays, neighbor);
                    if((int)(VECTOR(degrees)[neighbor]) <= k
                       && query_vertex_set.find(neighbor) != query_vertex_set.end()) {
                        turn_to_outer = true;
                        break;
                    }
                }
            }
            pending_vertices.set((size_t)vertex_idx);
            if(turn_to_outer)
                break;
            live_start = degree_starts[k];
        } else {
            bool checked = bfs_connected_queries(ori_graph, start_idx,
                                                 &query_vertex_set, &removed_vertices, &pending_vertices);
            if(checked) {
                removed_vertices |= pending_vertices;
                pending_vertices.reset();
            } else
                break;
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }
    auto vital_vertices = get_vital_component(ori_graph, start_idx, &removed_vertices);
    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

    if(vital_vertices.size() <= 20000) {
        vital_tool = induce_vertex_subgraph_mdc(ori_graph, &vital_vertices, induced_graph, induced_names, query_vertex_set);
        return 0;
    }
    vital_tool = vital_vertices;
    return 1;
}

unordered_map<int, double> compute_query_distance(igraph_t* graph, double & max_distance,
                                                  unordered_set<int>& new_query_set) {

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    unordered_map<int, double> result; //vertex-> query_distance
    boost::dynamic_bitset<> visited((size_t)igraph_vcount(graph));
    for(const auto & query_vertex: new_query_set) {
        queue<pair<int, int>*> bfs_queue;
        pair<int, int>* pair_element = new pair<int, int>(query_vertex, 0);
        bfs_queue.push(pair_element);

        visited.reset();
        visited.set((size_t)query_vertex);

        while(!bfs_queue.empty()) {
            auto bfs_e = bfs_queue.front();
            bfs_queue.pop();
            int source = bfs_e->first;
            int distance = bfs_e->second;
            delete bfs_e;

            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(!visited.test((size_t)neighbor)) {
                    visited.set((size_t)neighbor);
                    if(result.find(neighbor) == result.end()) {
                        result[neighbor] = (distance + 1) * (distance + 1);
                    } else {
                        result[neighbor] += (distance + 1) * (distance + 1);
                    }
                    if(max_distance < result[neighbor])
                        max_distance = result[neighbor];
                    pair_element = new pair<int, int>(neighbor, distance + 1);
                    bfs_queue.push(pair_element);
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return result;
}

unordered_map<int, double> compute_query_distance_copied(double & max_distance,
                                                         unordered_set<int> & new_query_set,
                                                         unordered_set<int> & vital_vertices,
                                                         igraph_vector_t & degrees,
                                                         boost::dynamic_bitset<>& ori_removed_vertices) {

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    igraph_degree(&copy_graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    for(int i = 0; i < igraph_vcount(&copy_graph); i++) {
        if(vital_vertices.find(i) == vital_vertices.end()) {
            ori_removed_vertices.set((size_t) i);
            VECTOR(degrees)[i] = 0;
            igraph_neighbors(&copy_graph, &neighbors, i, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];
                if(vital_vertices.find(neighbor) != vital_vertices.end())
                    VECTOR(degrees)[neighbor] -= 1;
            }
        }
    }

    unordered_map<int, double> result; //vertex-> query_distance
    boost::dynamic_bitset<> visited((size_t)igraph_vcount(&copy_graph));
    for(const auto & query_vertex: new_query_set) {
        queue<pair<int, int>*> bfs_queue;
        pair<int, int>* pair_element = new pair<int, int>(query_vertex, 0);
        bfs_queue.push(pair_element);

        visited.reset();
        visited.set((size_t)query_vertex);

        while(!bfs_queue.empty()) {
            auto bfs_e = bfs_queue.front();
            bfs_queue.pop();
            int source = bfs_e->first;
            int distance = bfs_e->second;
            delete bfs_e;

            igraph_neighbors(&copy_graph, &neighbors, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(!visited.test((size_t)neighbor) &&
                   !ori_removed_vertices.test((size_t)neighbor)) {
                    visited.set((size_t)neighbor);
                    if(result.find(neighbor) == result.end()) {
                        result[neighbor] = (distance + 1) * (distance + 1);
                    } else {
                        result[neighbor] += (distance + 1) * (distance + 1);
                    }
                    if(max_distance < result[neighbor])
                        max_distance = result[neighbor];
                    pair_element = new pair<int, int>(neighbor, distance + 1);
                    bfs_queue.push(pair_element);
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return result;
};

bool is_connected(igraph_t* graph,
                  int vertex_start_idx,
                  unordered_set<int>* vertex_idx_set,
                  boost::dynamic_bitset<>* removed_vertices) {
    if(removed_vertices->test((size_t)vertex_start_idx))
        return false;

    queue<int> bfs_queue;
    boost::dynamic_bitset<> visited(removed_vertices->size());
    int query_traveled_counter = 1;
    visited.set((size_t)vertex_start_idx);
    bfs_queue.push(vertex_start_idx);
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(!bfs_queue.empty()) {
        auto source = bfs_queue.front();
        bfs_queue.pop();
        igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            int neighbor = (int)VECTOR(neighbors)[i];
            if(!visited.test((size_t)neighbor)
               && !removed_vertices->test((size_t)neighbor)) {
                visited.set((size_t)neighbor);
                bfs_queue.push(neighbor);
                if(vertex_idx_set->find(neighbor) != vertex_idx_set->end())
                    query_traveled_counter++;
                if(query_traveled_counter >= vertex_idx_set->size()) {
                    igraph_vector_destroy(&neighbors);
                    return true;
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return false;
}

unordered_set<int> get_maximal_core_copy(igraph_vector_t & degrees,
                                         unordered_set<int>& query_set,
                                         boost::dynamic_bitset<>* elapsed_vertices) {
    int start_idx = *query_set.begin();

    auto degree_arrays = new int[igraph_vcount(&copy_graph)]();
    int max_degree = -1;
    for(int i = 0; i < igraph_vector_size(&degrees); i++) {
        auto degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree] += 1;
        if(degree > max_degree)
            max_degree = degree;
    }

    auto degree_starts = new int[max_degree + 1]();
    auto degree_ends = new int[max_degree + 1]();
    auto vertex_positions = new int[igraph_vcount(&copy_graph)]();

    for(int i = 1; i < max_degree + 1; i++) {
        degree_starts[i] = degree_starts[i - 1] + degree_arrays[i - 1];
        degree_ends[i] = degree_starts[i];
    }

    for(int i = 0; i < igraph_vcount(&copy_graph); i++) {
        auto vertex_degree = (int)(VECTOR(degrees)[i]);
        degree_arrays[degree_ends[vertex_degree]] = i;
        vertex_positions[i] = degree_ends[vertex_degree];
        degree_ends[vertex_degree]++;
    }

    int zero_end_idx = degree_ends[0];

    delete degree_ends;

    boost::dynamic_bitset<> pending_vertices((size_t)igraph_vcount(&copy_graph));
    boost::dynamic_bitset<> removed_vertices((size_t)igraph_vcount(&copy_graph));

    removed_vertices |= (*elapsed_vertices);

    auto k = (int)(VECTOR(degrees)[degree_arrays[zero_end_idx]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(&copy_graph);
    bool turn_to_outer = false;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    while(true) {
        if(live_start >= len_degree_arrays)
            break;
        int vertex_idx = degree_arrays[live_start];
        if((int)(VECTOR(degrees)[vertex_idx]) <= k) {
            del_degree_vertex(degree_starts, &degrees, vertex_idx);
            igraph_neighbors(&copy_graph, &neighbors, vertex_idx, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)(VECTOR(neighbors)[i]);
                if((int)(VECTOR(degrees)[neighbor]) > k) {
                    update_degree_vertex(degree_starts, &degrees, vertex_positions, degree_arrays, neighbor);
                    if((int)(VECTOR(degrees)[neighbor]) <= k
                       && query_set.find(neighbor) != query_set.end()) {
                        turn_to_outer = true;
                        break;
                    }
                }
            }
            pending_vertices.set((size_t)vertex_idx);
            if(turn_to_outer)
                break;
            live_start = degree_starts[k];
        } else {
            bool checked = bfs_connected_queries(&copy_graph, start_idx,
                                                 &query_set, &removed_vertices, &pending_vertices);
            if(checked) {
                removed_vertices |= pending_vertices;
                pending_vertices.reset();
            } else
                break;
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }
    auto vital_vertices = get_vital_component(&copy_graph, start_idx, &removed_vertices);
    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

    return vital_vertices;
}

//vertices orders in new graph
unordered_set<int> get_maximal_core(igraph_t* graph,
                                    unordered_set<int>& new_query_set,
                                    boost::dynamic_bitset<>* elapsed_vertices) {
    int start_idx = *new_query_set.begin();

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, igraph_vcount(graph));

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    for(int i = 0; i < igraph_vcount(graph); i++) {
        int degree = 0;
        if(!elapsed_vertices->test((size_t)i)) {
            igraph_neighbors(graph, &neighbors, i, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                int neighbor = (int)VECTOR(neighbors)[l];
                if(!elapsed_vertices->test((size_t)neighbor))
                    degree++;
            }
        }
        VECTOR(degrees)[i] = degree;
    }

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
    int zero_end_idx = degree_ends[0];
    delete degree_ends;

    boost::dynamic_bitset<> pending_vertices((size_t)igraph_vcount(graph));
    boost::dynamic_bitset<> removed_vertices((size_t)igraph_vcount(graph));

    removed_vertices |= (*elapsed_vertices);

    auto k = (int)(VECTOR(degrees)[degree_arrays[zero_end_idx]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(graph);
    bool turn_to_outer = false;

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
                    if((int)(VECTOR(degrees)[neighbor]) <= k
                       && new_query_set.find(neighbor) != new_query_set.end()) {
                        turn_to_outer = true;
                        break;
                    }
                }
            }
            pending_vertices.set((size_t)vertex_idx);
            if(turn_to_outer)
                break;
            live_start = degree_starts[k];
        } else {
            bool checked = bfs_connected_queries(graph, start_idx,
                                                 &new_query_set, &removed_vertices, &pending_vertices);
            if(checked) {
                removed_vertices |= pending_vertices;
                pending_vertices.reset();
            } else
                break;
            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }
    auto vital_vertices = get_vital_component(graph, start_idx, &removed_vertices);
    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);

    return vital_vertices;
}


// 0: live up to community size
// 1: still larger
// 2: not connected
int greedy_test(igraph_t * graph, vector<int>* names,
                double dist_limit, int upper_bound,
                unordered_set<int>& new_query_set,
                boost::dynamic_bitset<>* removed_vertices,
                unordered_map<int, double> & dist_info,
                unordered_set<int> & result) {
    removed_vertices->reset();
    for(const auto & dist_pair: dist_info) {
        if(dist_pair.second > dist_limit)
            removed_vertices->set((size_t)dist_pair.first);
    }
    if(!is_connected(graph,
                     *(new_query_set.begin()),
                     &new_query_set,
                     removed_vertices))
        return 2;
    auto temp_res = get_maximal_core(graph, new_query_set, removed_vertices);
    result.clear();
    for(const auto & element: temp_res) {
        auto name = names->at((unsigned long)element);
        result.insert(name);
    }

    if(temp_res.size() <= upper_bound)
        return 0;

    return 1;
}

int greedy_test_copy(double dist_limit, int upper_bound,
                     unordered_set<int>& query_set,
                     boost::dynamic_bitset<>* removed_vertices,
                     boost::dynamic_bitset<>* ori_removed_vertices,
                     unordered_map<int, double> & dist_info,
                     igraph_vector_t & new_degrees,
                     unordered_set<int> & result) {

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    removed_vertices->reset();
    for(const auto & dist_pair: dist_info) {
        if(dist_pair.second > dist_limit) {
            removed_vertices->set((size_t) dist_pair.first);
            VECTOR(new_degrees)[dist_pair.first] = 0;
            igraph_neighbors(&copy_graph, &neighbors, dist_pair.first, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];
                if(!ori_removed_vertices->test((size_t)neighbor)
                   && !removed_vertices->test((size_t)neighbor))
                    VECTOR(new_degrees)[neighbor] -= 1;
            }
        }
    }
    igraph_vector_destroy(&neighbors);

    *removed_vertices |= *ori_removed_vertices;
    if(!is_connected(&copy_graph,
                     *(query_set.begin()),
                     &query_set,
                     removed_vertices))
        return 2;

    result = get_maximal_core_copy(new_degrees, query_set, removed_vertices);

    if(result.size() <= upper_bound)
        return 0;

    return 1;
}

unordered_set<int> find_mdc(igraph_t * ori_graph,
                            igraph_vector_int_t * query_vertices, int upper_bound) {
    igraph_t induced_graph;
    vector<int> induced_names;

    unordered_set<int> vital_tool;

    int type = initial_decomp_max_core(query_vertices, ori_graph, &induced_graph, &induced_names, vital_tool);

    if(type == 0) {
        unordered_set<int> result;
        for(const auto & name:induced_names)
            result.insert(name);

        if(induced_names.size() <= upper_bound) {
            igraph_destroy(&induced_graph);
            return result;
        }

        double max_distance = -1.0;
        auto dist_info = compute_query_distance(&induced_graph, max_distance, vital_tool);
        double dist_lower = 0;
        double dist_upper = max_distance;
        double dist_limit;
        boost::dynamic_bitset<> removed_vertices((size_t)igraph_vcount(&induced_graph));

        dist_limit = (dist_lower + dist_upper) * 0.5;
        do {

            // 0: live up; 1: still larger; 2: not connected
            int mdc_errno = greedy_test(&induced_graph, &induced_names,
                                        dist_limit, upper_bound, vital_tool,
                                        &removed_vertices, dist_info, result);

            switch (mdc_errno) {
                case 0:
                    igraph_destroy(&induced_graph);
                    return result;
                case 1:
                    dist_upper = dist_limit;
                    break;
                case 2:
                    dist_lower = dist_limit;
                    break;
            }

            dist_limit = (dist_lower + dist_upper) * 0.5;
        } while (dist_limit > dist_lower && dist_limit < dist_upper);

        igraph_destroy(&induced_graph);
        return result;
    }

//    igraph_destroy(&induced_graph);

    unordered_set<int> result;
    unordered_set<int> new_query_set;
    for(int l = 0; l < igraph_vector_int_size(query_vertices); l++) {
        int query_vertex = VECTOR(*query_vertices)[l];
        new_query_set.insert(query_vertex);
    }

    igraph_vector_t degrees;
    igraph_vector_init(&degrees, 0);

    boost::dynamic_bitset<> ori_removed_vertices((size_t)igraph_vcount(&copy_graph));
    double max_distance = -1.0;
    auto dist_info = compute_query_distance_copied(max_distance,
                                                   new_query_set,
                                                   vital_tool, degrees,
                                                   ori_removed_vertices);

    double dist_lower = 0;
    double dist_upper = max_distance;
    double dist_limit;
    boost::dynamic_bitset<> removed_vertices((size_t)igraph_vcount(&copy_graph));

    dist_limit = (dist_lower + dist_upper) * 0.5;
    do {

        igraph_vector_t new_degrees;
        igraph_vector_init(&new_degrees, 0);
        igraph_vector_copy(&new_degrees, &degrees);

        // 0: live up; 1: still larger; 2: not connected
        int mdc_errno = greedy_test_copy(dist_limit, upper_bound, new_query_set,
                                         &removed_vertices, &ori_removed_vertices,
                                         dist_info, new_degrees, result);

        igraph_vector_destroy(&new_degrees);

        switch (mdc_errno) {
            case 0:
                igraph_vector_destroy(&degrees);
                return result;
            case 1:
                dist_upper = dist_limit;
                break;
            case 2:
                dist_lower = dist_limit;
                break;
        }

        dist_limit = (dist_lower + dist_upper) * 0.5;
    } while (dist_limit > dist_lower && dist_limit < dist_upper);


    igraph_vector_destroy(&degrees);
    return result;
}