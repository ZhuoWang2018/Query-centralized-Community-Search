//
// Created by wz on 18-2-7.
//

#include "steiner_tree.h"

double c = 3.0;
double alpha = 0.5;
double numerator;
igraph_t* steiner_graph;

int steiner_type;// 1-> ordinary; 2-> pruning; 3-> heuristic
int rank_limit = 50;

void init_steiner(igraph_t* g, int type) {
    steiner_graph = g;
    steiner_type = type;
    numerator = (igraph_ecount(g) + 0.0) / igraph_vcount(g) * 4;
}

void set_steiner_c(double data) {
    c = data;
}

void set_steiner_alpha(double data) {
    alpha = data;
}

unordered_map<int, unordered_set<int>*> initialize_connect_dict(igraph_vector_int_t * query_vertices) {
    unordered_map<int, unordered_set<int>*> connect_dict;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        auto key = VECTOR(*query_vertices)[i];
        auto value = new unordered_set<int>();
        value->insert((VECTOR(*query_vertices)[i]));
        connect_dict.insert({key, value});
    }
    return connect_dict;
}

double get_distance_bound(int query_node, bool query_connected,
                          unordered_map<int, double> & lookup_distance) {
    if(query_connected)
        return lookup_distance[query_node];
    return DBL_MAX;
}

double get_distance(igraph_vector_t* degrees) {
    auto minimum_degree = VECTOR(*degrees)[0];
    if(minimum_degree > VECTOR(*degrees)[1])
        minimum_degree = VECTOR(*degrees)[1];
    return (numerator + 0.0) / minimum_degree + 1;
}

bool add_connectivity(int source, int target, int query_vertices_number,
                      unordered_map<int, unordered_set<int>*>* connect_dict) {
    auto set1 = connect_dict->at(source);
    auto set2 = connect_dict->at(target);
    auto temp = set1;
    unordered_set<int> * p1;
    unordered_set<int> * p2;
    if(set1 != set2) {
        if(set1->size() > set2->size()) {
            p1 = set1;
            p2 = set2;
        } else {
            p1 = set2;
            p2 = set1;
        }
        for(const auto & element: *p2) {
            p1->insert(element);
            connect_dict->at(element) = p1;
        }
        delete p2;
        temp = p1;
    }
    return temp->size() >= query_vertices_number;
}

double dijkstra(igraph_t* graph, int origin, vector<double>* distances) {
    priority_queue<route*, vector<route*>, route_cmp> routes;

    unordered_set<int> visited;
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    auto vertices_size = igraph_vcount(graph);
    routes.push(new route(0, -1, origin, -1));
    while(!routes.empty()) {
        auto small_route = routes.top();
        routes.pop();

        auto distance = small_route->distance;
        auto newcomer = small_route->current_vertex_idx;

        delete small_route;
        if(visited.find(newcomer) != visited.end())
            continue;

        igraph_neighbors(graph, &neighbors, newcomer, IGRAPH_ALL);
        for(int k = 0; k < igraph_vector_size(&neighbors); k++) {
            auto neighbor = (int)VECTOR(neighbors)[k];
            if(visited.find(neighbor) == visited.end()) {
                igraph_integer_t eid;
                igraph_get_eid(graph, &eid, newcomer, neighbor, IGRAPH_UNDIRECTED, 0);
                auto new_distance = distances->at((unsigned long)eid) + distance;
                routes.push(new route(new_distance, -1, neighbor, -1));
            }
        }
        visited.insert(newcomer);
        if(visited.size() >= vertices_size) {
            igraph_vector_destroy(&neighbors);
            while(!routes.empty()) {
                auto route_to_del = routes.top();
                routes.pop();
                delete route_to_del;
            }
            return distance;
        }
    }

    igraph_vector_destroy(&neighbors);
    return DBL_MAX;
}



unordered_map<int, double> compute_distance_bounds(
        unordered_map<EDGE, route*, edge_hash, edge_equal>* pairs) {
    unordered_map<int, int> lookup_table;
    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);
    vector<double> distances;

    for(const auto & element: *pairs) {
        int source = get_induced_vertex(element.first->from, &lookup_table);
        int target = get_induced_vertex(element.first->to, &lookup_table);
        igraph_vector_push_back(&edges, source);
        igraph_vector_push_back(&edges, target);
        distances.push_back(element.second->distance);
    }

    igraph_t connected_graph;
    igraph_create(&connected_graph, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&edges);

    vector<int> names(lookup_table.size());
    for(const auto & element: lookup_table)
        names.at((unsigned long)element.second) = element.first;

    unordered_map<int, double> distance_bounds;
    for(int i = 0; i < names.size(); i++) {
        auto distance_bound = dijkstra(&connected_graph, i, &distances);
        distance_bounds[names[i]] = distance_bound;
    }

    igraph_destroy(&connected_graph);
    return distance_bounds;
};

void updateG1(route* from_route, route* to_route, igraph_vector_t* degrees,
              unordered_map<EDGE, route*, edge_hash, edge_equal>* triples,
              unordered_map<EDGE, route*, edge_hash, edge_equal>* pairs,
              int query_vertices_number, bool & query_connected,
              unordered_map<int, unordered_set<int>*>* connect_dict,
              unordered_map<int, double> & distance_bounds) {

    edge* surplus_pair = NULL;
    auto edge_pair = new edge(from_route->origin_vertex_idx,
                              to_route->origin_vertex_idx);

    if(triples->find(edge_pair) == triples->end())
        triples->insert({edge_pair, new route(DBL_MAX, -1, -1, -1)});
    else
        surplus_pair = edge_pair;

    auto distance = from_route->distance + to_route->distance
                    + get_distance(degrees);

    auto current_route = triples->at(edge_pair);
    if(current_route->distance > distance) {
        current_route->distance = distance;
        current_route->previous_vertex_idx = from_route->current_vertex_idx;
        current_route->current_vertex_idx = to_route->current_vertex_idx;
        if(!query_connected) {
            if(pairs->find(edge_pair) == pairs->end())
                pairs->insert({edge_pair, current_route});

            query_connected = add_connectivity(edge_pair->from, edge_pair->to,
                                               query_vertices_number, connect_dict);

            if(query_connected)
                distance_bounds = compute_distance_bounds(pairs);
        }
    }
    delete surplus_pair;
}

void updateG1(route* from_route, route* to_route, igraph_vector_t* degrees,
              unordered_map<EDGE, route*, edge_hash, edge_equal>* triples) {
    edge* surplus_pair = NULL;
    auto edge_pair = new edge(from_route->origin_vertex_idx,
                              to_route->origin_vertex_idx);
    if(triples->find(edge_pair) == triples->end())
        triples->insert({edge_pair, new route(DBL_MAX, -1, -1, -1)});
    else
        surplus_pair = edge_pair;

    auto distance = from_route->distance + to_route->distance
                    + get_distance(degrees);

    auto current_route = triples->at(edge_pair);
    if(current_route->distance > distance) {
        current_route->distance = distance;
        current_route->previous_vertex_idx = from_route->current_vertex_idx;
        current_route->current_vertex_idx = to_route->current_vertex_idx;
    }
    delete surplus_pair;
}

peak_node* get_top_element(unordered_map<int, handle_t>& dist_handles,
                           unordered_map<int, pres_handle_t>& pres_handles,
                           boost::heap::d_ary_heap<route*,
                                   boost::heap::mutable_<true>,
                                   boost::heap::arity<arity_number>,
                                   boost::heap::stable<true>,
                                   boost::heap::compare<route_cmp>>& routes,
                           boost::heap::d_ary_heap<pres_node*,
                                   boost::heap::mutable_<true>,
                                   boost::heap::arity<arity_number>,
                                   boost::heap::stable<true>,
                                   boost::heap::compare<pres_cmp>>& preses) {
    //vertex_idx -> peak_node
    unordered_map<int, peak_node*> top_elements;
    int rank = 1;

    int first_prestige_vertex = preses.top()->vertex_idx;

//    printf("routes size: %d\n",routes.size());
//    printf("prestige size: %d\n", preses.size());
    while(true) {
        if(!routes.empty()) {
            route* top_route_node = routes.top();
            int top_route_idx = top_route_node->current_vertex_idx;
            if(top_elements.find(top_route_idx) == top_elements.end()) {
                auto pn_dist = new peak_node();
                pn_dist->distance_rank = rank;
                pn_dist->distance = top_route_node->distance;
                top_elements[top_route_idx] = pn_dist;
                (*dist_handles[top_route_idx])->distance = DBL_MAX;
                routes.update(dist_handles[top_route_idx]);
            } else {
                auto pn_pres = top_elements[top_route_idx];

                pn_pres->distance_rank = rank;
                pn_pres->distance = top_route_node->distance;
                (*dist_handles[top_route_idx])->distance = -1.0;
                routes.update(dist_handles[top_route_idx]);

                (*pres_handles[top_route_idx])->prestige = DBL_MAX;
                preses.update(pres_handles[top_route_idx]);

                break;
            }


            pres_node* top_prestige_node = preses.top();
            int top_prestige_idx = top_prestige_node->vertex_idx;
            if(top_elements.find(top_prestige_idx) == top_elements.end()) {
                auto pn_pres = new peak_node();
                pn_pres->prestige_rank = rank;
                pn_pres->prestige = top_prestige_node->prestige;
                top_elements[top_prestige_idx] = pn_pres;
                (*pres_handles[top_prestige_idx])->prestige = -1.0;
                preses.update(pres_handles[top_prestige_idx]);
            } else {
                auto pn_dist = top_elements[top_prestige_idx];

                pn_dist->prestige_rank = rank;
                pn_dist->prestige = top_prestige_node->prestige;
                (*dist_handles[top_prestige_idx])->distance = -1.0;
                routes.update(dist_handles[top_prestige_idx]);

                (*pres_handles[top_prestige_idx])->prestige = DBL_MAX;
                preses.update(pres_handles[top_prestige_idx]);

                break;
            }
        }

        rank++;
        if(rank >= rank_limit) {
            auto pn_dist = top_elements[first_prestige_vertex];
            pn_dist->distance = (*dist_handles[first_prestige_vertex])->distance;
            pn_dist->distance_rank = rank_limit;
            (*dist_handles[first_prestige_vertex])->distance = -1.0;
            routes.update(dist_handles[first_prestige_vertex]);
            (*pres_handles[first_prestige_vertex])->prestige = DBL_MAX;
            preses.update(pres_handles[first_prestige_vertex]);
            break;
        }

    }

    peak_node * res = NULL;
    for(const auto & pair: top_elements) {
        int vertex_idx = pair.first;
        auto temp = pair.second;
        if(temp->prestige_rank > 0 && temp->distance_rank > 0) {
            res = temp;
        } else {
            if(temp->prestige_rank < 0) {//update dist_handles
                (*dist_handles[vertex_idx])->distance = temp->distance;
                routes.update(dist_handles[vertex_idx]);
            } else { // update pres_handles
                (*pres_handles[vertex_idx])->prestige = temp->prestige;
                preses.update(pres_handles[vertex_idx]);
            }
            delete temp;
        }
    }
    return res;
}

vector<route*> partition(const igraph_vector_int_t* query_vertices) {
    vector<route*> parlances(igraph_vcount(steiner_graph));

    boost::heap::d_ary_heap<route*,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<route_cmp>> routes;

    vector<handle_t> handles(igraph_vcount(steiner_graph));

    auto const null_handle = handle_t();

    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        int vertex = VECTOR(*query_vertices)[i];
        auto route_ele = new route(0, vertex, vertex, vertex);
        handles[vertex] = routes.push(route_ele);
    }

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    boost::dynamic_bitset<> visited((size_t)igraph_vcount(steiner_graph));

    igraph_vector_t degrees, seqs;
    igraph_vector_init(&degrees, 2);
    igraph_vector_init(&seqs, 2);

    while(!routes.empty()) {
        auto small_route = routes.top();
        routes.pop();

        auto newcomer = small_route->current_vertex_idx;
        auto distance = small_route->distance;
        auto ori_vertex_idx = small_route->origin_vertex_idx;

        parlances[newcomer] = small_route;
        igraph_neighbors(steiner_graph, &neighbors, newcomer, IGRAPH_ALL);

        for(int k = 0; k < igraph_vector_size(&neighbors); k++) {
            if(!visited.test((size_t)(VECTOR(neighbors)[k]))) {
                auto neighbor = (int)VECTOR(neighbors)[k];
                VECTOR(seqs)[0] = neighbor;
                VECTOR(seqs)[1] = newcomer;
                igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seqs),
                              IGRAPH_ALL, IGRAPH_LOOPS);
                auto new_distance = get_distance(&degrees) + distance;

                if(handles[neighbor] == null_handle) {
                    handles[neighbor] = routes.push(new route(new_distance, newcomer, neighbor, ori_vertex_idx));
                } else if(new_distance < (*(handles[neighbor]))->distance) {
                    route* temp_route = *(handles[neighbor]);
                    temp_route->distance = new_distance;
                    temp_route->previous_vertex_idx = newcomer;
                    temp_route->current_vertex_idx = neighbor;
                    temp_route->origin_vertex_idx = ori_vertex_idx;
                    routes.update(handles[neighbor]);
                }
            }
        }
        visited.set((size_t)newcomer);
    }

    igraph_vector_destroy(&neighbors);
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&seqs);
    return parlances;
}


vector<vector<int>*> construct_complete_graph_heuristic(
        igraph_vector_int_t* query_vertices,
        igraph_t* G1, igraph_vector_t* edge_distances) {

    int max_rank = -1;
    double min_prestige = 1.0;
    unordered_set<int> scanned;
    bool query_connected = false;
    unordered_map<int, int> lookup_table; // vertex_in_ori-> vertex_in_new
    unordered_map<int, double> distance_bounds; // query_vertex_in_ori-> query_distance
    unordered_map<int, pair<int, double>> lookup_relevance; // vertex_in_ori-> hyper_rank, relevance

    boost::heap::d_ary_heap<route*,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<route_cmp>> routes;
    unordered_map<int, handle_t> handles;

    boost::heap::d_ary_heap<pres_node*,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<pres_cmp>> preses;
    unordered_map<int, pres_handle_t> pres_handles;

    unordered_map<int, route*> parlances;
    unordered_map<EDGE, route*, edge_hash, edge_equal> triples;
    unordered_map<EDGE, route*, edge_hash, edge_equal> pairs;

    auto connect_dict = initialize_connect_dict(query_vertices);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    igraph_vector_t degrees, seq;
    igraph_vector_init(&degrees, 2);
    igraph_vector_init(&seq, 2);

    int query_vertices_number = igraph_vector_int_size(query_vertices);
    igraph_vector_t query_seqs;
    igraph_vector_init(&query_seqs, 0);
    for(int i = 0; i < query_vertices_number; i++)
        igraph_vector_push_back(&query_seqs, VECTOR(*query_vertices)[i]);

    igraph_vector_t query_degrees;
    igraph_vector_init(&query_degrees, query_vertices_number);
    igraph_degree(steiner_graph,
                  &query_degrees,
                  igraph_vss_vector(&query_seqs),
                  IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_vector_destroy(&query_seqs);

    for(int i = 0; i < query_vertices_number; i++) {
        auto vertex_idx = VECTOR(*query_vertices)[i];
        auto route_ele = new route(0, vertex_idx, vertex_idx, vertex_idx);
        handles[vertex_idx] = routes.push(route_ele);

        auto pres_ele = new pres_node(vertex_idx, 1.0, (int)VECTOR(query_degrees)[i]);
        pres_handles[vertex_idx] = preses.push(pres_ele);
    }
    igraph_vector_destroy(&query_degrees);

    while(!routes.empty()) {

        auto peak_vertex = get_top_element(handles, pres_handles, routes, preses);

        auto small_route = routes.top();
        routes.pop();

        auto small_prestige = preses.top();
        preses.pop();

        int temp_rank = peak_vertex->prestige_rank >= peak_vertex->distance_rank?
                        peak_vertex->prestige_rank: peak_vertex->distance_rank;

        small_route->distance = peak_vertex->distance;
        small_prestige->prestige = peak_vertex->prestige;

        delete peak_vertex;

        if(query_connected &&
           (temp_rank >= max_rank && small_prestige->prestige <= min_prestige)) {
            delete small_prestige;
            delete small_route;
            break;
        }

        auto newcomer = small_route->current_vertex_idx;
        auto distance = small_route->distance;
        auto ori_vertex_idx = small_route->origin_vertex_idx;

        parlances[newcomer] = small_route;

        if(!query_connected)
            lookup_relevance[newcomer] = make_pair(temp_rank, small_prestige->prestige);

        auto distance_bound = get_distance_bound(
                ori_vertex_idx, query_connected, distance_bounds);

        distance = distance * c;
        if(distance_bound <= distance) {
            scanned.insert(newcomer);
            delete small_prestige;
            continue;
        }

        igraph_neighbors(steiner_graph, &neighbors, newcomer, IGRAPH_ALL);
        auto new_prestige = (small_prestige->prestige + 0.0) / igraph_vector_size(&neighbors) * alpha;
        delete small_prestige;
        for(int k = 0; k < igraph_vector_size(&neighbors); k++) {
            auto neighbor = (int)VECTOR(neighbors)[k];
            if(scanned.find(neighbor) == scanned.end()) {
                VECTOR(seq)[0] = neighbor;
                VECTOR(seq)[1] = newcomer;
                igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seq),
                              IGRAPH_ALL, IGRAPH_LOOPS);
                auto new_distance = get_distance(&degrees) + distance;

                if(new_distance < distance_bound) {
                    if(handles.find(neighbor) == handles.end()) {
                        handles[neighbor] = routes.push(new route(new_distance, newcomer, neighbor, ori_vertex_idx));

                        igraph_degree(steiner_graph, &degrees, igraph_vss_1(neighbor), IGRAPH_ALL, IGRAPH_LOOPS);

                        pres_handles[neighbor] = preses.push(new pres_node(neighbor, new_prestige, VECTOR(degrees)[0]));
                    } else if(new_distance < (*(handles[neighbor]))->distance) {
                        route* temp_route = *(handles[neighbor]);
                        temp_route->distance = new_distance;
                        temp_route->previous_vertex_idx = newcomer;
                        temp_route->current_vertex_idx = neighbor;
                        temp_route->origin_vertex_idx = ori_vertex_idx;
                        routes.update(handles[neighbor]);

                        pres_node* temp_pres = *(pres_handles[neighbor]);
                        temp_pres->prestige += new_prestige;
                        preses.update(pres_handles[neighbor]);
                    }
                }
            } else {
                auto from_parlance = parlances.at(newcomer);
                auto to_parlance = parlances.at(neighbor);

                auto parent_from_idx = from_parlance->origin_vertex_idx;
                auto parent_to_idx = to_parlance->origin_vertex_idx;

                if(parent_from_idx == parent_to_idx)
                    continue;

                VECTOR(seq)[0] = neighbor;
                VECTOR(seq)[1] = newcomer;
                igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seq),
                              IGRAPH_ALL, IGRAPH_LOOPS);

                if(parent_from_idx < parent_to_idx)
                    updateG1(from_parlance, to_parlance,
                             &degrees, &triples, &pairs,
                             query_vertices_number, query_connected, &connect_dict, distance_bounds);
                else
                    updateG1(to_parlance, from_parlance,
                             &degrees, &triples, &pairs,
                             query_vertices_number, query_connected, &connect_dict, distance_bounds);
            }
        }
        scanned.insert(newcomer);
    }

    while(!routes.empty()) {
        auto small_route = routes.top();
        routes.pop();
        delete small_route;
        auto small_prestige = preses.top();
        preses.pop();
        delete small_prestige;
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&seq);

    igraph_vector_destroy(&neighbors);
    delete (*(connect_dict.begin())).second;

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    vector<vector<int>*> edges_paths;

    for(const auto & element: triples) {
        auto edge_source = get_induced_vertex(element.first->from, &lookup_table);
        auto edge_target = get_induced_vertex(element.first->to, &lookup_table);

        igraph_vector_push_back(&edges, edge_source);
        igraph_vector_push_back(&edges, edge_target);
        auto value = element.second;

        vector<int> * path = new vector<int>();
        vector<int> source_path;
        auto source_vertex_idx = value->previous_vertex_idx;
        auto parlance = parlances.at(source_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            source_path.push_back(parlance->current_vertex_idx);
            parlance = parlances.at(parlance->previous_vertex_idx);
        }
        source_path.push_back(parlance->current_vertex_idx);
        while(!source_path.empty()) {
            auto e = source_path.back();
            path->push_back(e);
            source_path.pop_back();
        }
        auto target_vertex_idx = value->current_vertex_idx;
        parlance = parlances.at(target_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            path->push_back(parlance->current_vertex_idx);
            parlance = parlances.at(parlance->previous_vertex_idx);
        }
        path->push_back(parlance->current_vertex_idx);

        edges_paths.push_back(path);
        igraph_vector_push_back(edge_distances, value->distance);
        delete element.first;
        delete element.second;
    }

    for(const auto & element: parlances)
        delete element.second;

    igraph_create(G1, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&edges);

    return edges_paths;
}

vector<vector<int>*> construct_complete_graph_prune(
        igraph_vector_int_t* query_vertices,
        igraph_t* G1, igraph_vector_t* edge_distances) {
    bool query_connected = false;
    unordered_map<int, int> lookup_table; // vertex_in_ori-> vertex_in_new
    unordered_map<int, double> distance_bounds; // query_vertex_in_ori-> query_distance

    unordered_set<int> scanned;
    boost::heap::d_ary_heap<route*,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<route_cmp>> routes;

    unordered_map<int, handle_t> handles;
    unordered_map<int, route*> parlances;
    unordered_map<EDGE, route*, edge_hash, edge_equal> triples;
    unordered_map<EDGE, route*, edge_hash, edge_equal> pairs;

    auto connect_dict = initialize_connect_dict(query_vertices);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    igraph_vector_t degrees, seq;
    igraph_vector_init(&degrees, 2);
    igraph_vector_init(&seq, 2);

    int query_vertices_number = igraph_vector_int_size(query_vertices);
    for(int i = 0; i < query_vertices_number; i++) {
        auto vertex_idx = VECTOR(*query_vertices)[i];
        auto route_ele = new route(0, vertex_idx, vertex_idx, vertex_idx);
        handles[vertex_idx] = routes.push(route_ele);
    }

    while(!routes.empty()) {
        auto small_route = routes.top();
        routes.pop();

        auto newcomer = small_route->current_vertex_idx;
        auto distance = small_route->distance;
        auto ori_vertex_idx = small_route->origin_vertex_idx;

        parlances[newcomer] = small_route;
        auto distance_bound = get_distance_bound(
                ori_vertex_idx, query_connected, distance_bounds);
        distance *= c;
        if(distance_bound <= distance) {
            scanned.insert(newcomer);
            continue;
        }

        igraph_neighbors(steiner_graph, &neighbors, newcomer, IGRAPH_ALL);
        for(int k = 0; k < igraph_vector_size(&neighbors); k++) {
            auto neighbor = (int)VECTOR(neighbors)[k];
            if(scanned.find(neighbor) == scanned.end()) {
                VECTOR(seq)[0] = neighbor;
                VECTOR(seq)[1] = newcomer;
                igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seq),
                              IGRAPH_ALL, IGRAPH_LOOPS);
                auto new_distance = get_distance(&degrees) + distance;

                if(new_distance < distance_bound) {
                    if(handles.find(neighbor) == handles.end()) {
                        handles[neighbor] = routes.push(
                                new route(new_distance, newcomer, neighbor, ori_vertex_idx));
                    } else if(new_distance < (*(handles[neighbor]))->distance) {
                        route* temp_route = *(handles[neighbor]);
                        temp_route->distance = new_distance;
                        temp_route->previous_vertex_idx = newcomer;
                        temp_route->current_vertex_idx = neighbor;
                        temp_route->origin_vertex_idx = ori_vertex_idx;
                        routes.update(handles[neighbor]);
                    }
                }

            } else {
                auto from_parlance = parlances.at(newcomer);
                auto to_parlance = parlances.at(neighbor);
                auto parent_from_idx = from_parlance->origin_vertex_idx;
                auto parent_to_idx = to_parlance->origin_vertex_idx;

                if(parent_from_idx == parent_to_idx)
                    continue;

                VECTOR(seq)[0] = neighbor;
                VECTOR(seq)[1] = newcomer;
                igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seq),
                              IGRAPH_ALL, IGRAPH_LOOPS);

                if(parent_from_idx < parent_to_idx)
                    updateG1(from_parlance, to_parlance, &degrees, &triples, &pairs,
                             query_vertices_number, query_connected, &connect_dict, distance_bounds);
                else
                    updateG1(to_parlance, from_parlance, &degrees, &triples, &pairs,
                             query_vertices_number, query_connected, &connect_dict, distance_bounds);

            }

        }
        scanned.insert(newcomer);
    }
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&seq);

    igraph_vector_destroy(&neighbors);
    delete (*(connect_dict.begin())).second;

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    vector<vector<int>*> edges_paths;
    for(const auto & element: triples) {
        auto edge_source = get_induced_vertex(element.first->from, &lookup_table);
        auto edge_target = get_induced_vertex(element.first->to, &lookup_table);

        igraph_vector_push_back(&edges, edge_source);
        igraph_vector_push_back(&edges, edge_target);
        auto value = element.second;

        vector<int> * path = new vector<int>();
        vector<int> source_path;
        auto source_vertex_idx = value->previous_vertex_idx;
        auto parlance = parlances.at(source_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            source_path.push_back(parlance->current_vertex_idx);
            parlance = parlances.at(parlance->previous_vertex_idx);
        }
        source_path.push_back(parlance->current_vertex_idx);
        while(!source_path.empty()) {
            auto e = source_path.back();
            path->push_back(e);
            source_path.pop_back();
        }
        auto target_vertex_idx = value->current_vertex_idx;
        parlance = parlances.at(target_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            path->push_back(parlance->current_vertex_idx);
            parlance = parlances.at(parlance->previous_vertex_idx);
        }
        path->push_back(parlance->current_vertex_idx);

        edges_paths.push_back(path);
        igraph_vector_push_back(edge_distances, value->distance);
        delete element.first;
        delete element.second;
    }

    for(const auto & element: parlances)
        delete element.second;

    igraph_create(G1, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&edges);

    return edges_paths;
}

vector<vector<int>*> construct_complete_graph_mehlhorn(
        igraph_vector_int_t* query_vertices,
        igraph_t* G1, igraph_vector_t* edge_distances) {
    bool query_connected = false;
    unordered_map<int, int> lookup_table;
    unordered_map<EDGE, route*, edge_hash, edge_equal> triples;

    auto parlances = partition(query_vertices);
    igraph_vector_t degrees, seqs;
    igraph_vector_init(&degrees, 2);
    igraph_vector_init(&seqs, 2);

    for(int i = 0; i < igraph_ecount(steiner_graph); i++) {
        auto from_vertex = IGRAPH_FROM(steiner_graph, i);
        auto to_vertex = IGRAPH_TO(steiner_graph, i);
        auto from_vertex_parlance = parlances[from_vertex];
        auto to_vertex_parlance = parlances[to_vertex];

        auto from_vertex_src = from_vertex_parlance->origin_vertex_idx;
        auto to_vertex_src = to_vertex_parlance->origin_vertex_idx;

        if(from_vertex_src == to_vertex_src)
            continue;

        VECTOR(seqs)[0] = from_vertex;
        VECTOR(seqs)[1] = to_vertex;
        igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seqs),
                      IGRAPH_ALL, IGRAPH_LOOPS);

        if(from_vertex < to_vertex)
            updateG1(parlances[from_vertex],
                     parlances[to_vertex], &degrees, &triples);
        else
            updateG1(parlances[to_vertex],
                     parlances[from_vertex], &degrees, &triples);
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&seqs);

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);
    vector<vector<int>*> edge_paths;

    for(const auto & element : triples) {
        auto edge_source = get_induced_vertex(element.first->from, &lookup_table);
        auto edge_target = get_induced_vertex(element.first->to, &lookup_table);
        igraph_vector_push_back(&edges, edge_source);
        igraph_vector_push_back(&edges, edge_target);
        auto value = element.second;

        vector<int> * path = new vector<int>();
        vector<int> source_path;
        auto source_vertex_idx = value->previous_vertex_idx;
        auto parlance = parlances.at((unsigned long)source_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            source_path.push_back(parlance->current_vertex_idx);
            parlance = parlances.at((unsigned long)parlance->previous_vertex_idx);
        }
        source_path.push_back(parlance->current_vertex_idx);
        while(!source_path.empty()) {
            auto e = source_path.back();
            path->push_back(e);
            source_path.pop_back();
        }
        auto target_vertex_idx = value->current_vertex_idx;
        parlance = parlances.at((unsigned long)target_vertex_idx);
        while(parlance->current_vertex_idx != parlance->origin_vertex_idx) {
            path->push_back(parlance->current_vertex_idx);
            parlance = parlances.at((unsigned long)parlance->previous_vertex_idx);
        }
        path->push_back(parlance->current_vertex_idx);

        edge_paths.push_back(path);
        igraph_vector_push_back(edge_distances, value->distance);
        delete element.first;
        delete element.second;
    }

    for(int i = 0; i < parlances.size(); i++) {
        auto element_route = parlances.at((unsigned long)i);
        delete element_route;
    }

    igraph_create(G1, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&edges);
    return edge_paths;
}

vector<int> replace_minimal_spanning_tree(igraph_vector_t * minimal_seqs, igraph_t * complete_graph,
                                          igraph_t * replaced_graph, vector<vector<int>*> * edge_paths,
                                          igraph_vector_t * replace_distances,
                                          igraph_vector_int_t* query_vertices,
                                          unordered_set<int> & query_nodes_in_replace) {
    igraph_vector_t replace_edges;
    igraph_vector_init(&replace_edges, 0);

    igraph_vector_t degrees, seq;
    igraph_vector_init(&degrees, 2);
    igraph_vector_init(&seq, 2);

    unordered_map<int, int> lookup_table;

    for(int i = 0; i < igraph_vector_size(minimal_seqs); i++) {
        auto edge_idx = (unsigned long)(VECTOR(*minimal_seqs)[i]);
        auto path = edge_paths->at(edge_idx);
        for(int k = 0; k < path->size() - 1; k++) {
            auto first_index = get_induced_vertex(path->at(k), &lookup_table);
            auto second_index = get_induced_vertex(path->at(k + 1), &lookup_table);
            igraph_vector_push_back(&replace_edges, first_index);
            igraph_vector_push_back(&replace_edges, second_index);

            VECTOR(seq)[0] = path->at((unsigned long)k);
            VECTOR(seq)[1] = path->at((unsigned long)(k + 1));
            igraph_degree(steiner_graph, &degrees, igraph_vss_vector(&seq),
                          IGRAPH_ALL, IGRAPH_LOOPS);

            auto edge_distance = get_distance(&degrees);
            igraph_vector_push_back(replace_distances, edge_distance);
        }
    }

    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&seq);

    igraph_create(replaced_graph, &replace_edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&replace_edges);
    vector<int> names((unsigned long)igraph_vcount(replaced_graph));
    for(const auto & element : lookup_table)
        names.at((unsigned long)element.second) = element.first;

    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        query_nodes_in_replace.insert(lookup_table[VECTOR(*query_vertices)[i]]);
    }

    return names;
}

bool visit(unordered_set<int>* tree_edge_set, unordered_set<int>* visited,
           unordered_map<int, unordered_set<int>*> *points,
           int steiner_start_point, igraph_t * replaced_graph, unordered_set<int> * steiner_indices) {
    visited->insert(steiner_start_point);

//    printf("elements in visit set\n");
//    for(const auto & element: *visited) {
//        printf("%d\n", element);
//    }
    bool has_steiner_point = false;
    auto edge_set = points->at(steiner_start_point);
    auto it = edge_set->begin();
    while(it != edge_set->end()) {
        auto edge_idx = *it;
        igraph_integer_t source;
        igraph_integer_t target;
        igraph_edge(replaced_graph, edge_idx, &source, &target);
        int neighbor;
        if(source == steiner_start_point)
            neighbor = target;
        else
            neighbor = source;

        if(visited->find(neighbor) != visited->end()) {
            it++;
            continue;
        }

        if(visit(tree_edge_set, visited, points, neighbor, replaced_graph, steiner_indices))
            has_steiner_point = true;
        else
            tree_edge_set->erase(edge_idx);
        it++;
    }
    if(steiner_indices->find(steiner_start_point) != steiner_indices->end())
        return true;
    return has_steiner_point;
}

unordered_set<int> get_steiner_edges(igraph_vector_t * minimal_seqs, igraph_t * replaced_graph,
                                     unordered_set<int> & query_nodes_in_replace, vector<int>* names) {
    unordered_map<int, unordered_set<int>*> points;
    unordered_set<int> tree_edge_set;
    for(int i = 0; i < igraph_vector_size(minimal_seqs); i++) {
        auto edge_idx = (int)VECTOR(*minimal_seqs)[i];
        tree_edge_set.insert(edge_idx);
    }
    unordered_set<int> visited;

    for(int i = 0; i < igraph_vector_size(minimal_seqs); i++) {
        auto edge_idx = (int)(VECTOR(*minimal_seqs)[i]);
        igraph_integer_t source;
        igraph_integer_t target;
        igraph_edge(replaced_graph, edge_idx, &source, &target);
        if(points.find(source) == points.end())
            points.insert({source, new unordered_set<int>()});
        if(points.find(target) == points.end())
            points.insert({target, new unordered_set<int>()});
        points.at(source)->insert(edge_idx);
        points.at(target)->insert(edge_idx);
    }

    auto steiner_start_point = *(query_nodes_in_replace.begin());
    visit(&tree_edge_set, &visited, &points,
          steiner_start_point, replaced_graph, &query_nodes_in_replace);

    unordered_set<int> edge_indices_in_ori_graph;
    auto it = tree_edge_set.begin();
    while(it != tree_edge_set.end()) {
        auto edge_idx = *it;
        igraph_integer_t source;
        igraph_integer_t target;
        igraph_edge(replaced_graph, edge_idx, &source, &target);
        auto from_in_ori_graph = names->at((unsigned long)source);
        auto to_in_ori_graph = names->at((unsigned long)target);
        igraph_integer_t eid;
        igraph_get_eid(steiner_graph, &eid, from_in_ori_graph, to_in_ori_graph, IGRAPH_UNDIRECTED, 0);
        edge_indices_in_ori_graph.insert(eid);
        it++;
    }

    for(const auto & element : points)
        delete element.second;

    return edge_indices_in_ori_graph;
}


unordered_set<int> get_steiner_tree(igraph_vector_int_t * query_vertices) {

    igraph_t G1;
    igraph_vector_t edge_distances;
    igraph_vector_init(&edge_distances, 0);

    vector<vector<int>*> edge_paths;
    if(steiner_type == 1) {
        edge_paths =
                construct_complete_graph_mehlhorn(query_vertices, &G1, &edge_distances);
    } else if (steiner_type == 2) {
        edge_paths =
                construct_complete_graph_prune(query_vertices, &G1, &edge_distances);
    } else if (steiner_type == 3) {
        edge_paths =
                construct_complete_graph_heuristic(query_vertices, &G1, &edge_distances);
    }

    igraph_vector_t steiner_edges;
    igraph_vector_init(&steiner_edges, 0);
    igraph_minimum_spanning_tree(&G1, &steiner_edges, &edge_distances);

    igraph_t replaced_graph;
    igraph_vector_t replace_distances;
    igraph_vector_init(&replace_distances, 0);

    unordered_set<int> query_nodes_in_replace;
    auto names = replace_minimal_spanning_tree(&steiner_edges, &G1, &replaced_graph,
                                               &edge_paths, &replace_distances,
                                               query_vertices, query_nodes_in_replace);
    igraph_destroy(&G1);
    igraph_vector_destroy(&steiner_edges);
    igraph_vector_destroy(&edge_distances);

    for(int i = 0; i < edge_paths.size(); i++)
        delete edge_paths.at(i);

    igraph_vector_t replace_steiner_edges;
    igraph_vector_init(&replace_steiner_edges, 0);
    igraph_minimum_spanning_tree(&replaced_graph, &replace_steiner_edges, &replace_distances);

    igraph_vector_destroy(&replace_distances);

    auto res = get_steiner_edges(
            &replace_steiner_edges, &replaced_graph, query_nodes_in_replace, &names);

    igraph_vector_destroy(&replace_steiner_edges);
    igraph_destroy(&replaced_graph);

    return res;
}
