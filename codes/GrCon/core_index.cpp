//
// Created by wz on 17-9-10.
//
#include "core_index.h"

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

void build_core_struct(igraph_t * graph, int core,
                       unordered_set<int> & pending_vertices,
                       unordered_set<int> & removed_vertices,
                       fstream & Hf) {

    printf("dealing with core:%d\n", core);

    igraph_vector_t neighbors;
    igraph_vector_init(&neighbors, 0);
    unordered_set<int> visited_pendings;

    int Horder = 0;
    for(const auto & element: pending_vertices) {
        if(visited_pendings.find(element) != visited_pendings.end())
            continue;

        queue<int> bfs_queue;
        bfs_queue.push(element);
        unordered_set<int> visited;
        visited.insert(element);
        visited_pendings.insert(element);

        while(!bfs_queue.empty()) {
            auto source = bfs_queue.front();
            bfs_queue.pop();

            igraph_neighbors(graph, &neighbors, source, IGRAPH_ALL);
            for(int l = 0; l < igraph_vector_size(&neighbors); l++) {
                auto neighbor = (int)VECTOR(neighbors)[l];
                if(removed_vertices.find(neighbor) == removed_vertices.end()
                   && visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                    if(pending_vertices.find(neighbor) != pending_vertices.end()
                            && visited_pendings.find(neighbor) == visited_pendings.end()) {
                        visited_pendings.insert(neighbor);
                    }
                }
            }
        }

        stringstream temp;
        temp<<core<<" "<<Horder<<"\n";
        for(const auto & e: visited) {
            temp<<e<<" ";
        }
        string temp_string = temp.str();
        temp_string = temp_string.substr(0, temp_string.size() - 1) + "\n";
        Hf<<temp_string;

        Horder++;
    }
    igraph_vector_destroy(&neighbors);
}

void core_decomposition(igraph_t * graph, fstream & Hf) {
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

    int zero_end_index = degree_ends[0];
    delete degree_ends;

    unordered_set<int> removed_vertices;
    unordered_set<int> pending_vertices;

    auto k = (int)(VECTOR(degrees)[degree_arrays[zero_end_index]]);
    int live_start = 0;
    int len_degree_arrays = igraph_vcount(graph);

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
                }
            }
            pending_vertices.insert(vertex_idx);
            live_start = degree_starts[k];
        } else {

            build_core_struct(graph, k, pending_vertices, removed_vertices, Hf);

            for(const auto & element: pending_vertices)
                removed_vertices.insert(element);
            pending_vertices.clear();

            k = (int)(VECTOR(degrees)[degree_arrays[live_start]]);
        }
    }

    igraph_vector_destroy(&neighbors);
    delete vertex_positions;
    delete degree_starts;
    delete degree_arrays;
    igraph_vector_destroy(&degrees);
}

