//
// Created by wz on 18-1-17.
//
#include "lctc.h"


igraph_t * g_lctc;
void init_lctc(igraph_t * g_data) {
    g_lctc = g_data;
    initialize_steiner_tree(g_lctc);
}

vector<int> expand_steiner_tree(unordered_set<int>* edge_indices, int upper_bound, int min_truss,
                                igraph_t * expanded_graph,
                                igraph_vector_int_t* query_vertices,
                                unordered_set<int> & query_nodes_in_expand) {

    unordered_map<int, int> lookup_table;

    unordered_set<int> edge_set;
    unordered_set<int> vertex_set;

    auto it_edge = edge_indices->begin();
    while(it_edge != edge_indices->end()) {
        auto edge_idx = *it_edge;
        edge_set.insert(edge_idx);
        int source = IGRAPH_FROM(g_lctc, edge_idx);
        int target = IGRAPH_TO(g_lctc, edge_idx);
        vertex_set.insert(source);
        vertex_set.insert(target);
        it_edge++;
    }

    queue<int> vertex_queue;
    auto it_vertex = vertex_set.begin();
    while(it_vertex != vertex_set.end()) {
        auto vertex_idx = *it_vertex;
        vertex_queue.push(vertex_idx);
        it_vertex++;
    }

    bool upper_bound_reached = false;
    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);

    while(vertex_queue.size() > 0 && !upper_bound_reached) {
        auto vertex_idx = vertex_queue.front();
        vertex_queue.pop();
        igraph_incident(g_lctc, &neighbor_eids, vertex_idx, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)VECTOR(neighbor_eids)[i];
            auto neighbor = IGRAPH_FROM(g_lctc, neighbor_eid);
            if(neighbor == vertex_idx)
                neighbor = IGRAPH_TO(g_lctc, neighbor_eid);
            if(edge_set.find(neighbor_eid) == edge_set.end()
               && (int)igraph_cattribute_EAN(g_lctc, "truss", neighbor_eid) >= min_truss) {
                edge_set.insert(neighbor_eid);
                if(vertex_set.find(neighbor) == vertex_set.end()) {
                    vertex_set.insert(neighbor);
                    if(vertex_set.size() >= upper_bound){
                        upper_bound_reached = true;
                        break;
                    }
                    vertex_queue.push(neighbor);
                }
            }
        }
    }
    igraph_vector_destroy(&neighbor_eids);

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    it_edge = edge_set.begin();
    while(it_edge != edge_set.end()) {
        auto ori_source = IGRAPH_FROM(g_lctc, *it_edge);
        auto ori_target = IGRAPH_TO(g_lctc, *it_edge);
        auto new_source = get_vertex_index(ori_source, lookup_table);
        auto new_target = get_vertex_index(ori_target, lookup_table);
        igraph_vector_push_back(&edges, new_source);
        igraph_vector_push_back(&edges, new_target);
        it_edge++;
    }

    igraph_create(expanded_graph, &edges, 0, IGRAPH_UNDIRECTED);
    vector<int> names((unsigned long)igraph_vcount(expanded_graph));
    for(const auto & element: lookup_table)
        names.at((unsigned long)element.second) = element.first;

    igraph_vector_destroy(&edges);

    for(int l = 0; l < igraph_vector_int_size(query_vertices); l++)
        query_nodes_in_expand.insert(lookup_table[VECTOR(*query_vertices)[l]]);

    return names;
}


vector<int> compute_graph_support(igraph_t * expanded_graph, unordered_map<int, unordered_map<int, int>> & neighbor_map) {

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    for(int i = 0; i < igraph_vcount(expanded_graph); i++) {
        igraph_incident(expanded_graph, &neighbor_eids, i, IGRAPH_ALL);
        for(int l = 0; l < igraph_vector_size(&neighbor_eids); l++) {
            auto neighbor_eid = (int)VECTOR(neighbor_eids)[l];
            int neighbor = IGRAPH_FROM(expanded_graph, neighbor_eid);
            if(neighbor == i)
                neighbor = IGRAPH_TO(expanded_graph, neighbor_eid);
            neighbor_map[i][neighbor] = neighbor_eid;
        }
    }
    igraph_vector_destroy(&neighbor_eids);
    vector<int> supports((unsigned long)igraph_ecount(expanded_graph));

    for(int i = 0; i < supports.size(); i++) {
        if(i % 100000 == 0 && i > 0)
            printf("i am still computing supports, on %d\n", i);
        int source = IGRAPH_FROM(expanded_graph, i);
        int target = IGRAPH_TO(expanded_graph, i);

        int support = 0;

        if(neighbor_map[source].size() < neighbor_map[target].size()) {
            for(const auto & element: neighbor_map[source]) {
                if(neighbor_map[target].find(element.first)
                   != neighbor_map[target].end())
                    support++;
            }
        } else {
            for(const auto & element: neighbor_map[target]) {
                if(neighbor_map[source].find(element.first)
                   != neighbor_map[source].end())
                    support++;
            }
        }
        supports[i] = support;
    }

    return supports;
}

void del_support_index(int * support_starts, vector<int>* supports, int edge_idx) {
    auto support = supports->at((unsigned long)edge_idx);
    support_starts[support] ++;
    supports->at((unsigned long)edge_idx)--;
}

void update_support_index(int * support_starts, vector<int>* supports,
                          int * edge_positions, int * support_arrays, int edge_idx, int k) {
    if(supports->at((unsigned long)edge_idx) > k) {
        auto support = supports->at((unsigned long)edge_idx);
        auto position_1 = support_starts[support];
        auto position_2 = edge_positions[edge_idx];
        if(position_1 != position_2) {
            auto temp = support_arrays[position_1];
            support_arrays[position_1] = support_arrays[position_2];
            support_arrays[position_2] = temp;
            edge_positions[edge_idx] = position_1;
            edge_positions[support_arrays[position_2]] = position_2;
        }
        support_starts[support] += 1;
        supports->at((unsigned long)edge_idx) --;
    }
}

bool bfs_connected_queries(igraph_t * graph, int vertex_start_idx, unordered_set<int>* vertex_idx_set,
                           unordered_set<int>* removed_edges, unordered_set<int>* pending_edges) {
    queue<int> bfs_queue;
    unordered_set<int> visited;
    int query_traveled_counter = 1;
    visited.insert(vertex_start_idx);
    bfs_queue.push(vertex_start_idx);
    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    while(!bfs_queue.empty()) {
        auto source = bfs_queue.front();
        bfs_queue.pop();
        igraph_incident(graph, &neighbor_eids, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)VECTOR(neighbor_eids)[i];
            int neighbor = IGRAPH_FROM(graph, neighbor_eid);
            if(neighbor == source)
                neighbor = IGRAPH_TO(graph, neighbor_eid);
            if(visited.find(neighbor) == visited.end()
               && removed_edges->find(neighbor_eid) == removed_edges->end()
               && pending_edges->find(neighbor_eid) == pending_edges->end()) {
                visited.insert(neighbor);
                bfs_queue.push(neighbor);
                if(vertex_idx_set->find(neighbor) != vertex_idx_set->end())
                    query_traveled_counter ++;
                if(query_traveled_counter >= vertex_idx_set->size()) {
                    igraph_vector_destroy(&neighbor_eids);
                    return true;
                }
            }
        }
    }
    igraph_vector_destroy(&neighbor_eids);
    return false;
}

unordered_set<int> get_vital_component(igraph_t * graph, int vertex_start_idx, unordered_set<int> * removed_edges) {
    queue<int> vital_queue;
    unordered_set<int> visited;
    visited.insert(vertex_start_idx);
    vital_queue.push(vertex_start_idx);

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    while(!vital_queue.empty()) {
        auto source = vital_queue.front();
        vital_queue.pop();

        igraph_incident(graph, &neighbor_eids, source, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
            auto neighbor = IGRAPH_FROM(graph, neighbor_eid);
            if(neighbor == source)
                neighbor = IGRAPH_TO(graph, neighbor_eid);
            if(visited.find(neighbor) == visited.end()
               && removed_edges->find(neighbor_eid) == removed_edges->end()) {
                visited.insert(neighbor);
                vital_queue.push(neighbor);
            }
        }
    }
    igraph_vector_destroy(&neighbor_eids);
    return visited;
}


void induce_vertex_subgraph(igraph_t * old_graph, vector<int>* names,
                            igraph_t * induced_graph, vector<int>* induced_names,
                            unordered_set<int> * vital_vertices) {
    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);

    unordered_map<int, int> induced_map;
    unordered_set<int> visited_edges;

    for(const auto & element: *vital_vertices){
        igraph_incident(old_graph, &neighbor_eids, element, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
            auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
            if(visited_edges.find(neighbor_eid) == visited_edges.end()) {
                int from = IGRAPH_FROM(old_graph, neighbor_eid);
                int to = IGRAPH_TO(old_graph, neighbor_eid);
                if(vital_vertices->find(from) != vital_vertices->end()
                   && vital_vertices->find(to) != vital_vertices->end()) {
                    auto from_vertex_idx = get_vertex_index(from, induced_map);
                    auto to_vertex_idx = get_vertex_index(to, induced_map);
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

    induced_names->resize(vital_vertices->size());
    for(const auto &element: induced_map) {
        induced_names->at((unsigned long)element.second) = names->at((unsigned long)element.first);
    }
}

vector<int> compute_maximal_truss_graph(igraph_t * expanded_graph, igraph_vector_int_t* query_vertices, vector<int>* names,
                                        igraph_t * maximal_truss_graph, int * maximal_truss,
                                        vector<int>* preserved_names,
                                        unordered_map<int, unordered_map<int, int>> & upper_neighbor_map,
                                        unordered_set<int> & query_nodes_in_expand) {

    int bfs_start_vertex_idx = *(query_nodes_in_expand.begin());

    unordered_map<int, unordered_map<int, int>> neighbor_map;
//    compute support of each edge
    vector<int> supports = compute_graph_support(expanded_graph, neighbor_map);

//    sort supports in the bucket way
    auto support_arrays = new int[igraph_ecount(expanded_graph)]();

    int max_support = -1;
    for(int i = 0; i < igraph_ecount(expanded_graph); i++) {
        auto support = supports.at((unsigned long)i);
        support_arrays[support] += 1;
        if(support > max_support)
            max_support = support;
    }

    auto support_starts = new int[max_support + 1]();
    auto support_ends = new int[max_support + 1]();
    auto edge_positions = new int[igraph_ecount(expanded_graph)]();

    for(int i = 1; i < max_support + 1; i++) {
        support_starts[i] = support_starts[i - 1] + support_arrays[i - 1];
        support_ends[i] = support_starts[i];
    }

    for(int i = 0; i < igraph_ecount(expanded_graph); i++) {
        auto edge_support = supports.at((unsigned long)i);
        support_arrays[support_ends[edge_support]] = i;
        edge_positions[i] = support_ends[edge_support];
        support_ends[edge_support]++;
    }
    delete support_ends;

    unordered_set<int> removed_edges;
    unordered_set<int> pending_edges;

    auto k = supports.at((unsigned long)support_arrays[0]);
    int live_start = 0;
    int len_support_arrays = igraph_ecount(expanded_graph);

    while(true){
        if(live_start >= len_support_arrays)
            break;
        int edge_idx = support_arrays[live_start];
        if(supports.at((unsigned long)edge_idx) <= k) {
            del_support_index(support_starts, &supports, edge_idx);
            auto source = IGRAPH_FROM(expanded_graph, edge_idx);
            auto target = IGRAPH_TO(expanded_graph, edge_idx);

            int edge_index1, edge_index2;
            int smaller, bigger;
            if(neighbor_map[source].size() < neighbor_map[target].size()) {
                smaller = source;
                bigger = target;
            } else {
                smaller = target;
                bigger = source;
            }

            for(const auto & element: neighbor_map[smaller]) {
                if(neighbor_map[bigger].find(element.first)
                   != neighbor_map[bigger].end()) {
                    edge_index1 = element.second;
                    edge_index2 = neighbor_map[bigger][element.first];
                    if(removed_edges.find(edge_index1) == removed_edges.end()
                       && removed_edges.find(edge_index2) == removed_edges.end()
                       && pending_edges.find(edge_index1) == pending_edges.end()
                       && pending_edges.find(edge_index2) == pending_edges.end()) {
                        update_support_index(support_starts, &supports, edge_positions,
                                             support_arrays, edge_index1, k);
                        update_support_index(support_starts, &supports, edge_positions,
                                             support_arrays, edge_index2, k);
                    }
                }
            }

            pending_edges.insert(edge_idx);
            live_start = support_starts[k];
        } else {
            auto checked = bfs_connected_queries(expanded_graph, bfs_start_vertex_idx,
                                                 &query_nodes_in_expand, &removed_edges, &pending_edges);
            if(checked) {
                auto it = pending_edges.begin();
                while(it != pending_edges.end()) {
                    removed_edges.insert(*it);
                    it++;
                }
                pending_edges.clear();
            } else
                break;
            k = supports.at((unsigned long)support_arrays[live_start]);
        }
    }
    auto vital_vertices = get_vital_component(expanded_graph, bfs_start_vertex_idx, &removed_edges);

    delete edge_positions;
    delete support_starts;
    delete support_arrays;


    induce_vertex_subgraph(expanded_graph, names, maximal_truss_graph, preserved_names, &vital_vertices);

    auto maximal_truss_supports = compute_graph_support(maximal_truss_graph, upper_neighbor_map);

    *maximal_truss = k;

    return maximal_truss_supports;
}

int compute_max_distance_lctc(igraph_t* graph, unordered_set<int> & query_set, unordered_set<int> * removed_edges) {
    int max_distance = -1;
    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    for(const auto & element: query_set) {
        queue<pair<int, int>*> bfs_queue;
        pair<int, int>* pair_element = new pair<int, int>(element, 0);
        bfs_queue.push(pair_element);

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            auto bfs_e = bfs_queue.front();
            bfs_queue.pop();
            int source = bfs_e->first;
            int distance = bfs_e->second;
            delete bfs_e;
            igraph_incident(graph, &neighbor_eids, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
                auto neighbor_eid = (int)(VECTOR(neighbor_eids)[i]);
                if(removed_edges->find(neighbor_eid) == removed_edges->end()) {
                    auto neighbor = IGRAPH_FROM(graph, neighbor_eid);
                    if(neighbor == source)
                        neighbor = IGRAPH_TO(graph, neighbor_eid);
                    if(visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        if(distance + 1 > max_distance)
                            max_distance = distance + 1;
                        pair_element = new pair<int, int>(neighbor, distance + 1);
                        bfs_queue.push(pair_element);
                    }
                }
            }
        }
    }
    igraph_vector_destroy(&neighbor_eids);
    return max_distance;
}

void filter_edge(igraph_t * graph, int vertex_idx, unordered_set<int>* removed_edges,
                 unordered_set<int>* pending_edges, igraph_vector_t * eids) {
    igraph_incident(graph, eids, vertex_idx, IGRAPH_ALL);
    for(int i = 0; i < igraph_vector_size(eids); i++) {
        int eid = (int)(VECTOR(*eids)[i]);
        if(removed_edges->find(eid) == removed_edges->end())
            pending_edges->insert(eid);
    }
}

unordered_set<int> compute_removed_edges(igraph_t * graph, unordered_set<int> * vertex_idx_set, unordered_set<int> * removed_edges, int * extreme_distance) {
    unordered_map<int, pair<int, int>*> dist_dict;

    auto it_vertex = vertex_idx_set->begin();

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    while(it_vertex != vertex_idx_set->end()) {
        queue<pair<int, int>*> bfs_queue;
        pair<int, int>* pair_element = new pair<int, int>(*it_vertex, 0);
        bfs_queue.push(pair_element);

        unordered_set<int> visited;
        visited.insert(*it_vertex);

        while(!bfs_queue.empty()) {
            auto bfs_e = bfs_queue.front();
            bfs_queue.pop();
            int source = bfs_e->first;
            int distance = bfs_e->second;
            delete bfs_e;
            igraph_incident(graph, &neighbor_eids, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbor_eids); i++) {
                auto neighbor_eid = (int)VECTOR(neighbor_eids)[i];
                int neighbor = IGRAPH_FROM(graph, neighbor_eid);
                if(neighbor == source)
                    neighbor = IGRAPH_TO(graph, neighbor_eid);
                if(visited.find(neighbor) == visited.end()) {
                    if(removed_edges->find(neighbor_eid) == removed_edges->end()) {
                        visited.insert(neighbor);
                        if(dist_dict.find(neighbor) == dist_dict.end()) {
                            auto dist_dict_element = new pair<int, int>(0, 0);
                            dist_dict.insert({neighbor, dist_dict_element});
                        }
                        auto dist_info = dist_dict.at(neighbor);
                        if(distance + 1 > dist_info->first)
                            dist_info->first = distance + 1;
                        dist_info->second += (distance + 1);
                        pair_element = new pair<int, int>(neighbor, distance + 1);
                        bfs_queue.push(pair_element);
                    }
                }
            }
        }
        it_vertex++;
    }
    igraph_vector_destroy(&neighbor_eids);

    igraph_vector_t eids;
    igraph_vector_init(&eids, 0);
    unordered_set<int> pending_edges;
    int max_distance = -1;
    int max_journey = -1;
    auto dict_it = dist_dict.begin();
    while(dict_it != dist_dict.end()) {
        auto distance = (*dict_it).second->first;
        if(max_distance < distance) {
            max_distance = distance;
            pending_edges.clear();
            max_journey = (*dict_it).second->second;
            filter_edge(graph, (*dict_it).first, removed_edges, &pending_edges, &eids);
        } else if(max_distance == distance) {
            if(max_journey < (*dict_it).second->second){
                pending_edges.clear();
                max_journey = (*dict_it).second->second;
                filter_edge(graph, (*dict_it).first, removed_edges, &pending_edges, &eids);
            } else if(max_journey == (*dict_it).second->second) {
                filter_edge(graph, (*dict_it).first, removed_edges, &pending_edges, &eids);
            }
        }
        delete (*dict_it).second;
        dict_it++;
    }
    igraph_vector_destroy(&eids);
    *extreme_distance = max_distance;
    return pending_edges;
}

unordered_set<int> maintain_maximal_truss(igraph_t * graph, unordered_set<int>* edges_to_del,
                                          unordered_set<int>* removed_edges, vector<int>* supports,
                                          int maximal_truss, unordered_map<int, unordered_map<int, int>> & neighbor_map) {
    unordered_set<int> deleted_edges;
    while(edges_to_del->size() > 0) {
        int edge_idx = *(edges_to_del->begin());
        edges_to_del->erase(edge_idx);
        deleted_edges.insert(edge_idx);
        int source = IGRAPH_FROM(graph, edge_idx);
        int target = IGRAPH_TO(graph, edge_idx);

        int edge_index1, edge_index2;
        int smaller, bigger;

        if(neighbor_map[source].size() < neighbor_map[target].size()) {
            smaller = source;
            bigger = target;
        } else {
            smaller = target;
            bigger = source;
        }

        for(const auto & element: neighbor_map[smaller]) {
            if(neighbor_map[bigger].find(element.first)
               != neighbor_map[bigger].end()) {
                edge_index1 = element.second;
                edge_index2 = neighbor_map[bigger][element.first];
                if(removed_edges->find(edge_index1) == removed_edges->end()
                   && removed_edges->find(edge_index2) == removed_edges->end()
                   && deleted_edges.find(edge_index1) == deleted_edges.end()
                   && deleted_edges.find(edge_index2) == deleted_edges.end()) {
                    supports->at((unsigned long) edge_index1) -= 1;
                    supports->at((unsigned long) edge_index2) -= 1;
                    if (supports->at((unsigned long) edge_index1) < maximal_truss)
                        edges_to_del->insert(edge_index1);
                    if (supports->at((unsigned long) edge_index2) < maximal_truss)
                        edges_to_del->insert(edge_index2);
                }
            }
        }
    }
    return deleted_edges;
}

unordered_set<int> bulk_deletion(igraph_t * graph, igraph_vector_int_t * query_vertices, int maximal_truss,
                                 vector<int>* supports, vector<int>* preserved_names,
                                 unordered_map<int, unordered_map<int, int>> & neighbor_map) {
    vector<int> vertex_idx_array;
    unordered_set<int> vertex_idx_set;

    unordered_set<int> query_vertex_idx_set;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        int query_vertex_idx = VECTOR(*query_vertices)[i];
        query_vertex_idx_set.insert(query_vertex_idx);
    }

    for(int i = 0; i < preserved_names->size(); i++) {
        int vertex_name = preserved_names->at((unsigned long)i);
        if(query_vertex_idx_set.find(vertex_name) != query_vertex_idx_set.end())
            vertex_idx_array.push_back(i);
    }

    for(int i = 0; i < vertex_idx_array.size(); i++) {
        int vertex_idx = vertex_idx_array.at((unsigned long)i);
        vertex_idx_set.insert(vertex_idx);
    }

    queue<unordered_set<int>*> removed_queue;
    unordered_set<int> * removed_edges = new unordered_set<int>();
    auto last_distance = compute_max_distance_lctc(graph, vertex_idx_set, removed_edges);
    unordered_set<int> * last_removed_edges = removed_edges;
    int no_use;

    while(true){
        auto pending_edges = compute_removed_edges(graph, &vertex_idx_set, removed_edges, &no_use);
        pending_edges = maintain_maximal_truss(graph, &pending_edges, removed_edges, supports, maximal_truss, neighbor_map);

        auto connected = bfs_connected_queries(graph, vertex_idx_array.at(0), &vertex_idx_set, removed_edges, &pending_edges);
        if(!connected)
            break;
        unordered_set<int>* new_removed_edges = new unordered_set<int>();

        auto it_edge = pending_edges.begin();
        while(it_edge != pending_edges.end()) {
            new_removed_edges->insert(*it_edge);
            it_edge++;
        }
        it_edge = removed_edges->begin();
        while(it_edge != removed_edges->end()) {
            new_removed_edges->insert(*it_edge);
            it_edge++;
        }

        removed_queue.push(removed_edges);
        removed_edges = new_removed_edges;

        int distance = compute_max_distance_lctc(graph, vertex_idx_set, removed_edges);

        if(connected && last_distance >= distance) {
            last_distance = distance;
            last_removed_edges = new_removed_edges;
            while(!removed_queue.empty()) {
                auto element = removed_queue.front();
                removed_queue.pop();
                delete element;
            }
        }
    }

    unordered_set<int> community_points;

    auto rest_vertex_indices = get_vital_component(graph, vertex_idx_array.at(0), last_removed_edges);

    auto it_vertex = rest_vertex_indices.begin();
    while(it_vertex != rest_vertex_indices.end()) {
        int name = preserved_names->at((unsigned long)(*it_vertex));
        community_points.insert(name);
        it_vertex++;
    }

    delete removed_edges;
    while(!removed_queue.empty()) {
        auto element = removed_queue.front();
        removed_queue.pop();
        delete element;
    }

    return community_points;

}

unordered_set<int> fetch_community(
        igraph_vector_int_t * query_vertices, int upper_bound) {

    unordered_set<int> steiner_edge_set;
    if(igraph_vector_int_size(query_vertices) == 1) {
        int query_vertex = VECTOR(*query_vertices)[0];
        igraph_vector_t eids;
        igraph_vector_init(&eids, 0);
        igraph_incident(g_lctc, &eids, query_vertex, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&eids); i++) {
            int eid = (int)VECTOR(eids)[i];
            steiner_edge_set.insert(eid);
        }
        igraph_vector_destroy(&eids);
    } else {
        steiner_edge_set = get_steiner_tree(query_vertices);
    }

    auto min_truss = INT32_MAX;

    auto it = steiner_edge_set.begin();
    while(it != steiner_edge_set.end()) {
        auto edge_truss = (int)igraph_cattribute_EAN(g_lctc, "truss", *it);
        if(min_truss > edge_truss)
            min_truss = edge_truss;
        it++;
    }

    igraph_t expanded_graph;
    unordered_set<int> query_nodes_in_expand;
    vector<int> names = expand_steiner_tree(&steiner_edge_set,
                                            upper_bound, min_truss, &expanded_graph,
                                            query_vertices, query_nodes_in_expand);


    igraph_t maximal_truss_graph;
    int maximal_truss;

    vector<int> preserved_names;
    unordered_map<int, unordered_map<int, int>> neighbor_map;
    auto maximal_truss_supports = compute_maximal_truss_graph(&expanded_graph, query_vertices, &names, &maximal_truss_graph,
                                                              &maximal_truss, &preserved_names, neighbor_map, query_nodes_in_expand);
    igraph_destroy(&expanded_graph);

    unordered_set<int> community_points = bulk_deletion(&maximal_truss_graph, query_vertices,
                                                        maximal_truss, &maximal_truss_supports,
                                                        &preserved_names, neighbor_map);

    igraph_destroy(&maximal_truss_graph);

    return community_points;
}
