//
// Created by wz on 17-9-10.
//

#include "online.h"
#include "graph_creator.h"

unordered_map<int, int> vertex_indices;//vertex-> core_idx
unordered_map<int, vector<boost::dynamic_bitset<>>> core_indices;//{core-> List<Horders>}


void split_ego(string & str, vector<int> & result) {
    string::size_type pos;
    result.clear();
    string pattern = " ";
    str += pattern;
    unsigned long size = (int)str.size();

    for(unsigned long i = 0; i < size; i++) {
        pos = str.find(" ", i);
        if(pos < size) {
            string s = str.substr(i, pos - i);
            int number;
            sscanf(s.data(), "%d", &number);
            result.push_back(number);
            i = pos + pattern.size() - 1;
        }
    }
}

void split_igraph(string str, igraph_vector_int_t * result) {
    string::size_type pos;
    igraph_vector_int_clear(result);
    string pattern = " ";
    str += pattern;
    unsigned long size = (int)str.size();

    for(unsigned long i = 0; i < size; i++) {
        pos = str.find(" ", i);
        if(pos < size) {
            string s = str.substr(i, pos - i);
            int number;
            sscanf(s.data(), "%d", &number);
            igraph_vector_int_push_back(result, number);
            i = pos + pattern.size() - 1;
        }
    }
}


void load_index(fstream & index_f, igraph_t * graph) {
    string s;
    vector<int> temp;
    while(getline(index_f, s)) {
        split_ego(s, temp);
        int core = temp[0];

        getline(index_f, s);
        split_ego(s, temp);
        boost::dynamic_bitset<> Hindex((size_t)igraph_vcount(graph));
        for(const auto & vertex: temp) {
            vertex_indices[vertex] = core;
            Hindex.set((size_t)vertex);
        }

        core_indices[core].push_back(Hindex);
    }
}

//0: not contained, test the next connected component
//1: not contained, just skip the current core
//2: contained, just use this
int test_connected(boost::dynamic_bitset<> & Hindex,
                   igraph_vector_int_t * query_vertices) {
    int query_counter = 0;
    for(int i = 0; i < igraph_vector_int_size(query_vertices); i++) {
        auto query_vertex = VECTOR(*query_vertices)[i];
        if(Hindex.test((size_t)query_vertex))
            query_counter++;
        else if(query_counter > 0)
            return 1;
        else
            return 0;
    }
    return 2;
}

boost::dynamic_bitset<> retrieval(igraph_t* graph, igraph_vector_int_t * query_vertices, int & k_core) {
    int min_core = vertex_indices[VECTOR(*query_vertices)[0]];
    for(int l = 1; l < igraph_vector_int_size(query_vertices); l++) {
        int query_vertex = VECTOR(*query_vertices)[l];
        if(min_core > vertex_indices[query_vertex])
            min_core = vertex_indices[query_vertex];
    }

    while(true) {
        for(int l = 0; l < core_indices[min_core].size(); l++) {
            boost::dynamic_bitset<> bitset = core_indices[min_core][l];
            int error_type = test_connected(bitset, query_vertices);
            if(error_type == 0)
                continue;
            else if(error_type == 1)
                break;
            else {
                k_core = min_core;
                return bitset;
            }
        }
        min_core--;
    }
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

unordered_set<int> get_neighbors(boost::dynamic_bitset<> & smaller_graph, int vertex, igraph_t * graph) {
    unordered_set<int> filter_neighbors;

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    igraph_neighbors(graph, &neighbors, vertex, IGRAPH_ALL);
    for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
        auto neighbor = (int)VECTOR(neighbors)[l];
        if(smaller_graph.test((size_t)neighbor))
            filter_neighbors.insert(neighbor);
    }
    igraph_vector_destroy(&neighbors);
    return filter_neighbors;
}


unordered_set<int> get_grcon_community(igraph_t* graph,
                                       igraph_vector_int_t* query_vertices) {
    int k_core;
    boost::dynamic_bitset<> Hcandidate = retrieval(graph, query_vertices, k_core);

//    unordered_set<int> smaller_graph_nodes;
//    int start_vertex_idx = (int)Hcandidate.find_first();
//    while(start_vertex_idx != Hcandidate.npos) {
//        smaller_graph_nodes.insert(start_vertex_idx);
//        start_vertex_idx = Hcandidate.find_next(start_vertex_idx);
//    }

    unordered_set<int> ori_query_set;
    for(int l = 0; l < igraph_vector_int_size(query_vertices); l++)
        ori_query_set.insert(VECTOR(*query_vertices)[l]);

    // build a steiner tree from the smaller graph

//    igraph_t induced_graph;
//    vector<int> induced_names;
//
//    auto new_query_set = induce_vertex_subgraph_mdc(graph, &smaller_graph_nodes, &induced_graph, &induced_names, ori_query_set);

    unordered_set<int> steiner_points;

    if(igraph_vector_int_size(query_vertices) > 1) {
        initialize_steiner_tree(graph, &Hcandidate);

        auto steiner_edges = get_steiner_tree(query_vertices);


        for(const auto & steiner_edge: steiner_edges) {
            auto point = IGRAPH_FROM(graph, steiner_edge);
            steiner_points.insert(point);

            point = IGRAPH_TO(graph, steiner_edge);
            steiner_points.insert(point);
        }
    } else {
        steiner_points = ori_query_set;
    }

//    printf("steiner is done\n");


//    igraph_destroy(&induced_graph);

    unordered_map<int, int> node_scores;
    unordered_map<int, unordered_set<int>> neighborCache;

    boost::heap::d_ary_heap<NodeRank,
            boost::heap::mutable_<true>,
            boost::heap::arity<arity_number>,
            boost::heap::stable<true>,
            boost::heap::compare<NodeRank_cmp>> queue;
    unordered_map<int, handle_score_t> handles;

    unordered_set<int> newQueueNodes;
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    for(const auto & steiner_point: steiner_points) {
        int degree = 0;
        igraph_neighbors(graph, &neighbors, steiner_point, IGRAPH_ALL);
        for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
            auto neighbor = (int)VECTOR(neighbors)[l];
            if(steiner_points.find(neighbor) != steiner_points.end())
                degree++;
        }
        int score = k_core - degree;
        if(score < 0)
            score = 0;
        node_scores[steiner_point] = score;

        if(neighborCache.find(steiner_point)
                == neighborCache.end())
            neighborCache[steiner_point] = get_neighbors(Hcandidate, steiner_point, graph);

        for(const auto & neighbor: neighborCache[steiner_point]) {
            if(steiner_points.find(neighbor)
                    == steiner_points.end())
                newQueueNodes.insert(neighbor);
        }
    }
    igraph_vector_destroy(&neighbors);

    bool hasAllScoresZeroFlag = true;
    for (const auto & pair : node_scores) {
        if(pair.second > 0) {
            hasAllScoresZeroFlag = false;
            break;
        }
    }

    if(hasAllScoresZeroFlag) {
        unordered_set<int> result;
        for(const auto & pair: node_scores)
            result.insert(pair.first);
    }

    for(const auto & newNode: newQueueNodes) {
        if(neighborCache.find(newNode)
                == neighborCache.end())
            neighborCache[newNode] = get_neighbors(Hcandidate, newNode, graph);

        NodeRank queueNode = NodeRank(newNode, k_core, node_scores, neighborCache[newNode]);
        handles[newNode] = queue.push(queueNode);
    }

    while(!queue.empty()) {
        NodeRank nodeRank = queue.top();
        queue.pop();
        handles.erase(nodeRank.node);

        int score = k_core - nodeRank.solutionDegree;
        if(score < 0)
            score = 0;
        node_scores[nodeRank.node] = score;

        if(neighborCache.find(nodeRank.node)
                == neighborCache.end())
            neighborCache[nodeRank.node] = get_neighbors(Hcandidate, nodeRank.node, graph);

        auto neighbors = neighborCache[nodeRank.node];

        for(const auto & neighbor: neighbors) {
            if(node_scores.find(neighbor)
                    != node_scores.end()) {
                if(node_scores[neighbor] > 0) {
                    int new_score = node_scores[neighbor] - 1;
                    if(new_score < 0)
                        new_score = 0;
                    node_scores[neighbor] = new_score;

                    if(new_score == 0) {
                        if(neighborCache.find(neighbor)
                                == neighborCache.end())
                            neighborCache[neighbor] = get_neighbors(Hcandidate, neighbor, graph);

                        auto neighborNeighbors = neighborCache[neighbor];

                        for(const auto & neighborNeighbor: neighborNeighbors) {
                            if(handles.find(neighborNeighbor)
                                    != handles.end()) {
                                (*handles[neighborNeighbor]).updateRankForSaturation();
                                queue.update(handles[neighborNeighbor]);
                            }
                        }
                    }
                }
            } else {

                if(handles.find(neighbor) != handles.end()) {
                    (*handles[neighbor]).updateRanksForAddition(nodeRank.node, k_core, node_scores);
                    queue.update(handles[neighbor]);
                } else {
                    if(neighborCache.find(neighbor)
                            == neighborCache.end())
                        neighborCache[neighbor] = get_neighbors(Hcandidate, neighbor, graph);
                    auto neighborsNeighbors = neighborCache[neighbor];
                    NodeRank queueNode(neighbor, k_core, node_scores, neighborsNeighbors);
                    handles[neighbor] = queue.push(queueNode);
                }
            }
        }

        for (const auto & pair : node_scores) {
            if (pair.second == 0) {
                hasAllScoresZeroFlag = true;
            } else {
                hasAllScoresZeroFlag = false;
                break;
            }
        }

        if (hasAllScoresZeroFlag) {
            break;
        }

        if(node_scores.size() > 100)
            break;
    }

    unordered_set<int> result;
    for(const auto & pair: node_scores)
        result.insert(pair.first);
    return result;
}


void grcon_test(string & file_name,
                string & index_name,
                string & query_file,
                string & result_file) {
    igraph_t g;

    fstream query_f(query_file, ios::in);
    fstream result_f(result_file, ios::out);

    get_simple_graph(file_name.data(), &g);

    fstream index_f(index_name, ios::in);
    load_index(index_f, &g);
    index_f.close();

    int community_counter = 0;

    string s;
    igraph_vector_int_t query_points;
    igraph_vector_int_init(&query_points, 0);
    unordered_set<int> fetched_points;

    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed = start - start;
    while(getline(query_f, s)) {
        split_igraph(s, &query_points);

        start = chrono::system_clock::now();
        fetched_points = get_grcon_community(&g, &query_points);
        end = chrono::system_clock::now();
        printf("one query has been answered!\n");
        elapsed += (end - start);

        stringstream s_stream;
        for(const auto & element: fetched_points) {
            s_stream<< element<< " ";
        }
        string result = s_stream.str();
        result_f<<result.substr(0, result.size() - 1)<<"\n";

        community_counter++;
        if(community_counter % 10 == 0) {
            printf("i am still working, on %d\n", community_counter);
            printf("time elapsed: %fs\n", elapsed.count());
        }
    }
    igraph_vector_int_destroy(&query_points);

    query_f.close();
    result_f.close();
    igraph_destroy(&g);
}
