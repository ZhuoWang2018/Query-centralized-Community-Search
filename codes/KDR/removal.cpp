//
// Created by wz on 18-2-7.
//

#include "removal.h"

unordered_set<int> maintain_k_core(igraph_t* g,
                                   unordered_set<int>* vertices_to_del,
                                   unordered_set<int>* removed_vertices,
                                   igraph_vector_t* degrees, int k_core) {
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    unordered_set<int> deleted_vertices;
    while(vertices_to_del->size() > 0) {
        int vertex_idx = *(vertices_to_del->begin());
        vertices_to_del->erase(vertex_idx);
        deleted_vertices.insert(vertex_idx);
        igraph_neighbors(g, &neighbors, vertex_idx, IGRAPH_ALL);
        for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
            auto neighbor = (int)(VECTOR(neighbors)[i]);
            if(removed_vertices->find(neighbor) == removed_vertices->end()
               && deleted_vertices.find(neighbor) == deleted_vertices.end()) {
                VECTOR(*degrees)[neighbor]--;
                if((int)(VECTOR(*degrees)[neighbor]) < k_core)
                    vertices_to_del->insert(neighbor);
            }
        }
    }
    igraph_vector_destroy(&neighbors);
    return deleted_vertices;
}




unordered_set<int> get_max_query_dist_vertices(igraph_t* g,
                                               unordered_set<int>* query_nodes,
                                               unordered_set<int>* removed_vertices,
                                               int & max_distance, bool bulk) {
    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);

    unordered_map<int, int> vertices_info; // vertex-> max_query_distance
    for(const auto & element: *query_nodes) {
        queue<pair<int,int>> bfs_queue;
        bfs_queue.push(make_pair(element, 0));

        unordered_set<int> visited;
        visited.insert(element);

        while(!bfs_queue.empty()) {
            int source = bfs_queue.front().first;
            int distance = bfs_queue.front().second;
            bfs_queue.pop();

            igraph_neighbors(g, &neighbors, source, IGRAPH_ALL);
            for(int i = 0; i < igraph_vector_size(&neighbors); i++) {
                auto neighbor = (int)VECTOR(neighbors)[i];
                if(removed_vertices->find(neighbor) == removed_vertices->end()
                   && visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    if(vertices_info.find(neighbor) == vertices_info.end()
                       || distance + 1 > vertices_info[neighbor])
                        vertices_info[neighbor] = distance + 1;
                    bfs_queue.push(make_pair(neighbor, distance + 1));
                }
            }
        }
    }
    igraph_vector_destroy(&neighbors);

    unordered_set<int> vertices;
    max_distance = -1;
    for(const auto & element: vertices_info) {
        auto distance = element.second;
        if(max_distance < distance) {
            max_distance = distance;
            vertices.clear();
            vertices.insert(element.first);
        } else if(max_distance == distance)
            vertices.insert(element.first);
    }
    if(bulk)
        return vertices;
    unordered_set<int> single_vertices;
    single_vertices.insert((*vertices.begin()));
    return single_vertices;
}