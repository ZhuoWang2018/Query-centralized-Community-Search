//
// Created by wz on 17-9-12.
//
#include "grcon_steiner_tree.h"

igraph_t * g;
boost::dynamic_bitset<>* H;
unordered_map<int, int> point_names;

double MAX_DISTANCE = 1.79769e+308;

void initialize_steiner_tree(igraph_t* gtemp,
                             boost::dynamic_bitset<> * Htemp) {
    g = gtemp;
    H = Htemp;
}

vector<route*> dijkstra_binary_global_bfs(const igraph_vector_int_t * vertex_idx_array, unordered_set<int> & visited_edges) {

    vector<route*> parlances((unsigned long)igraph_vcount(g));

    boost::heap::d_ary_heap<route*,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<route_cmp>> routes;

    vector<handle_t> handles((unsigned long)igraph_vcount(g));
    auto const null_handle = handle_t();

    for(int i = 0; i < igraph_vector_int_size(vertex_idx_array); i++) {
        int vertex_idx = VECTOR(*vertex_idx_array)[i];
        auto route_ele = new route(0, 0, vertex_idx, vertex_idx, vertex_idx);
        auto l = routes.push(route_ele);
        handles[vertex_idx] = l;
    }

    igraph_vector_t neighbor_eids;
    igraph_vector_init(&neighbor_eids, 0);
    boost::dynamic_bitset<> visited((size_t)igraph_vcount(g));

//    unordered_set<int> handles_put;

    while(!routes.empty()) {

        auto small_route = routes.top();

        routes.pop();

        auto newcomer = small_route->current_vertex_idx;
        auto distance = small_route->distance;
        auto ori_vertex_idx = small_route->origin_vertex_idx;

        parlances[(unsigned long)newcomer] = small_route;
        igraph_incident(g, &neighbor_eids, newcomer, IGRAPH_ALL);

        for(int k = 0; k < igraph_vector_size(&neighbor_eids); k++) {
            auto neighbor_eid = (int)VECTOR(neighbor_eids)[k];
            int neighbor = IGRAPH_FROM(g, neighbor_eid);
            if(neighbor == newcomer)
                neighbor = IGRAPH_TO(g, neighbor_eid);
            if(!visited.test((size_t)neighbor) && H->test((size_t)neighbor)) {

                visited_edges.insert(neighbor_eid);

                auto new_distance = distance + 1;

                if(handles[neighbor] == null_handle) {
                    handles[neighbor] = routes.push(new route(0, new_distance, newcomer, neighbor, ori_vertex_idx));
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
    igraph_vector_destroy(&neighbor_eids);
    return parlances;
}

void fill_triple_dict_bfs(route* from_vertex_route, route* to_vertex_route, int current_min_truss,
                          unordered_map<EDGE, route*, edge_hash, edge_equal>* triple_dict) {
    edge * surplus_pair = NULL;
    auto edge_pair = new edge(from_vertex_route->origin_vertex_idx,
                              to_vertex_route->origin_vertex_idx);

    auto distance = from_vertex_route->distance + to_vertex_route->distance + 1;

    if(triple_dict->find(edge_pair) == triple_dict->end()) {
        triple_dict->insert({edge_pair, new route(0, MAX_DISTANCE, -1, -1, -1)});
    } else {
        surplus_pair = edge_pair;
    }

    auto current_route = triple_dict->at(edge_pair);
    if(current_route->distance > distance) {
        current_route->distance = distance;
        current_route->min_truss = 0;
        current_route->previous_vertex_idx = from_vertex_route->current_vertex_idx;
        current_route->current_vertex_idx = to_vertex_route->current_vertex_idx;
    }
    delete surplus_pair;
}

int get_vertex_index(int vertex_name) {
    if(point_names.find(vertex_name) == point_names.end()) {
        auto dict_size = point_names.size();
        point_names.insert({vertex_name, dict_size});
        return (int)dict_size;
    }
    return point_names.at(vertex_name);
}

vector<vector<int>*> construct_complete_graph_global_bfs(igraph_vector_int_t * vertex_idx_array, igraph_t * new_graph, igraph_vector_t * edge_distances) {
    point_names.clear();
    unordered_map<EDGE, route*, edge_hash, edge_equal>triple_dict;


    unordered_set<int> visited_edges;
    auto parlances = dijkstra_binary_global_bfs(vertex_idx_array, visited_edges);

//    for(int i = 0; i < igraph_ecount(g); i++) {
    for(const auto & i : visited_edges) {

        auto from_vertex_idx = IGRAPH_FROM(g, i);
        auto to_vertex_idx = IGRAPH_TO(g, i);
        auto from_vertex_parlance = parlances.at((unsigned long)from_vertex_idx);
        auto to_vertex_parlance = parlances.at((unsigned long)to_vertex_idx);

        if(from_vertex_parlance == NULL || to_vertex_parlance == NULL)
            continue;
        auto from_vertex_src = from_vertex_parlance->origin_vertex_idx;
        auto to_vertex_src = to_vertex_parlance->origin_vertex_idx;

        if(from_vertex_src == to_vertex_src)
            continue;

        if(from_vertex_src < to_vertex_src) {
            fill_triple_dict_bfs(parlances.at((unsigned long)from_vertex_idx),
                                 parlances.at((unsigned long)to_vertex_idx),
                                 0,
                                 &triple_dict);
        } else {
            fill_triple_dict_bfs(parlances.at((unsigned long)to_vertex_idx),
                                 parlances.at((unsigned long)from_vertex_idx),
                                 0,
                                 &triple_dict);
        }
    }

    igraph_vector_t edges;
    igraph_vector_init(&edges, 0);

    vector<vector<int>*> edge_paths;

    for(const auto & element : triple_dict) {
        auto edge_source = get_vertex_index(element.first->from);
        auto edge_target = get_vertex_index(element.first->to);
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


    igraph_create(new_graph, &edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&edges);
    point_names.clear();
    return edge_paths;
}

vector<int> replace_minimal_spanning_tree_bfs(igraph_vector_t * minimal_seqs, igraph_t * complete_graph,
                                              igraph_t * replaced_graph, vector<vector<int>*> * edge_paths,
                                              igraph_vector_t * replace_distances) {
    point_names.clear();
    igraph_vector_t replace_edges;
    igraph_vector_init(&replace_edges, 0);

    for(int i = 0; i < igraph_vector_size(minimal_seqs); i++) {
        auto edge_idx = (unsigned long)(VECTOR(*minimal_seqs)[i]);
        auto path = edge_paths->at(edge_idx);
        for(int k = 0; k < path->size() - 1; k++) {
            auto first_index = get_vertex_index(path->at((unsigned long) k));
            auto second_index = get_vertex_index(path->at((unsigned long) (k + 1)));
            igraph_vector_push_back(&replace_edges, first_index);
            igraph_vector_push_back(&replace_edges, second_index);
            igraph_integer_t eid;
            igraph_get_eid(g, &eid, path->at((unsigned long) k), path->at((unsigned long) (k + 1)),
                           IGRAPH_UNDIRECTED, 0);
            auto edge_distance = 1;
            igraph_vector_push_back(replace_distances, edge_distance);
        }
    }
    igraph_create(replaced_graph, &replace_edges, 0, IGRAPH_UNDIRECTED);
    igraph_vector_destroy(&replace_edges);
    vector<int> names((unsigned long)igraph_vcount(replaced_graph));
    for(const auto & element : point_names)
        names.at((unsigned long)element.second) = element.first;

    return names;
}

bool visit(unordered_set<int>* tree_edge_set, unordered_set<int>* visited,
           unordered_map<int, unordered_set<int>*> *points,
           int steiner_start_point, igraph_t * replaced_graph, unordered_set<int> * steiner_indices) {
    visited->insert(steiner_start_point);

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
                                     igraph_vector_int_t * vertex_idx_array, vector<int>* names) {
    unordered_map<int, unordered_set<int>*> points;
    unordered_set<int> steiner_indices;
    unordered_set<int> tree_edge_set;
    for(int i = 0; i < igraph_vector_size(minimal_seqs); i++) {
        auto edge_idx = (int)VECTOR(*minimal_seqs)[i];
        tree_edge_set.insert(edge_idx);
    }
    unordered_set<int> visited;
    for(int i = 0; i < igraph_vector_int_size(vertex_idx_array); i++) {
        auto vertex_idx = point_names.at(VECTOR(*vertex_idx_array)[i]);
        steiner_indices.insert(vertex_idx);
    }

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

    auto steiner_start_point = point_names.at(VECTOR(*vertex_idx_array)[0]);
    visit(&tree_edge_set, &visited, &points, steiner_start_point, replaced_graph, &steiner_indices);

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
        igraph_get_eid(g, &eid, from_in_ori_graph, to_in_ori_graph, IGRAPH_UNDIRECTED, 0);
        edge_indices_in_ori_graph.insert(eid);
        it++;
    }
    point_names.clear();

    for(const auto & element : points)
        delete element.second;

    return edge_indices_in_ori_graph;
}

unordered_set<int> get_steiner_tree(igraph_vector_int_t * vertex_idx_array) {
    igraph_t new_graph;
    igraph_vector_t edge_distances;
    igraph_vector_init(&edge_distances, 0);

    vector<vector<int>*> edge_paths;
    edge_paths = construct_complete_graph_global_bfs(vertex_idx_array, &new_graph, &edge_distances);

    igraph_vector_t steiner_edges;
    igraph_vector_init(&steiner_edges, 0);
    igraph_minimum_spanning_tree(&new_graph, &steiner_edges, &edge_distances);

    igraph_t replaced_graph;
    igraph_vector_t replace_distances;
    igraph_vector_init(&replace_distances, 0);
    vector<int> names;

    names = replace_minimal_spanning_tree_bfs(&steiner_edges, &new_graph, &replaced_graph,
                                                  &edge_paths, &replace_distances);
    igraph_destroy(&new_graph);
    igraph_vector_destroy(&steiner_edges);
    igraph_vector_destroy(&edge_distances);

    for(int i = 0; i < edge_paths.size(); i++)
        delete edge_paths.at((unsigned long)i);

    igraph_vector_t replace_steiner_edges;
    igraph_vector_init(&replace_steiner_edges, 0);
    igraph_minimum_spanning_tree(&replaced_graph, &replace_steiner_edges, &replace_distances);

    igraph_vector_destroy(&replace_distances);

    auto res = get_steiner_edges(&replace_steiner_edges, &replaced_graph, vertex_idx_array, &names);

    igraph_vector_destroy(&replace_steiner_edges);
    igraph_destroy(&replaced_graph);

    return res;
}