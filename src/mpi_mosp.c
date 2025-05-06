#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>
#include <metis.h>

// Structure to represent an edge
typedef struct {
    int u;
    int v;
    int weight;
} Edge;

// Structure for adjacency list node
typedef struct AdjListNode {
    int dest;
    int weight;
    struct AdjListNode* next;
} AdjListNode;

// Structure for adjacency list
typedef struct {
    AdjListNode* head;
} AdjList;

// Structure for graph
typedef struct {
    int num_vertices;
    int num_edges;
    AdjList* array;
} Graph;

// Structure for SOSP tree
typedef struct {
    int* parent;
    int* distance;
    bool* marked;
} SOSPTree;

// Structure for graph partition
typedef struct {
    int start_vertex;
    int end_vertex;
    int* global_to_local;
    int* local_to_global;
    int* boundary_vertices;
    int num_boundary;
} Partition;

// Function prototypes
AdjListNode* newAdjListNode(int dest, int weight);
Graph* createGraph(int n);
void addEdge(Graph* graph, int src, int dest, int weight);
SOSPTree* initSOSPTree(int n, int source);
void freeSOSPTree(SOSPTree* tree);
void freeGraph(Graph* graph);
void SOSP_Update_Parallel(Graph* graph, SOSPTree* tree, Edge* inserted_edges, int num_inserted, int rank, int size, Partition* partition);
void computeInitialSOSP_Parallel(Graph* graph, SOSPTree* tree, int source, int rank, int size, Partition* partition);
Graph* readGraph(const char* filename);
Edge* createSmartRandomEdges(Graph* graph, int source, int target, int num_edges);
void MOSP_Update_2obj_Parallel(Graph* graph1, Graph* graph2, SOSPTree* tree1, SOSPTree* tree2, 
                              Edge* inserted_edges1, Edge* inserted_edges2, 
                              int num_inserted1, int num_inserted2, int source,
                              int rank, int size, Partition* partition);
Partition* partitionGraphMETIS(Graph* graph, int rank, int size);
void freePartition(Partition* partition);

// Implementation of existing functions (unchanged)
AdjListNode* newAdjListNode(int dest, int weight) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    if (newNode == NULL) {
        fprintf(stderr, "Memory allocation failed for adjacency list node\n");
        exit(EXIT_FAILURE);
    }
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = NULL;
    return newNode;
}

Graph* createGraph(int n) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    if (graph == NULL) {
        fprintf(stderr, "Memory allocation failed for graph\n");
        exit(EXIT_FAILURE);
    }
    
    graph->num_vertices = n;
    graph->num_edges = 0;
    
    graph->array = (AdjList*)malloc(n * sizeof(AdjList));
    if (graph->array == NULL) {
        fprintf(stderr, "Memory allocation failed for adjacency lists\n");
        free(graph);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        graph->array[i].head = NULL;
    }
    
    return graph;
}

void addEdge(Graph* graph, int src, int dest, int weight) {
    AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
    
    graph->num_edges++;
}

SOSPTree* initSOSPTree(int n, int source) {
    SOSPTree* tree = (SOSPTree*)malloc(sizeof(SOSPTree));
    if (tree == NULL) {
        fprintf(stderr, "Memory allocation failed for SOSP tree\n");
        exit(EXIT_FAILURE);
    }
    
    tree->parent = (int*)malloc(n * sizeof(int));
    tree->distance = (int*)malloc(n * sizeof(int));
    tree->marked = (bool*)malloc(n * sizeof(bool));
    
    if (tree->parent == NULL || tree->distance == NULL || tree->marked == NULL) {
        fprintf(stderr, "Memory allocation failed for SOSP tree arrays\n");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        tree->distance[i] = INT_MAX;
        tree->parent[i] = -1;
        tree->marked[i] = false;
    }
    
    tree->distance[source] = 0;
    
    return tree;
}

void freeSOSPTree(SOSPTree* tree) {
    if (tree != NULL) {
        free(tree->parent);
        free(tree->distance);
        free(tree->marked);
        free(tree);
    }
}

void freeGraph(Graph* graph) {
    if (graph != NULL) {
        for (int i = 0; i < graph->num_vertices; i++) {
            AdjListNode* current = graph->array[i].head;
            while (current) {
                AdjListNode* temp = current;
                current = current->next;
                free(temp);
            }
        }
        free(graph->array);
        free(graph);
    }
}

// Graph partitioning using METIS
Partition* partitionGraphMETIS(Graph* graph, int rank, int size) {
    int n = graph->num_vertices;
    Partition* partition = (Partition*)malloc(sizeof(Partition));
    if (partition == NULL) {
        fprintf(stderr, "Memory allocation failed for partition\n");
        exit(EXIT_FAILURE);
    }
    
    // Only process 0 does the partitioning
    idx_t nVertices = n;
    idx_t nConstraints = 1;
    idx_t nParts = size;
    idx_t objval;
    
    // Arrays needed for METIS
    idx_t* xadj = NULL;
    idx_t* adjncy = NULL;
    idx_t* part = NULL;
    
    if (rank == 0) {
        // Convert graph to CSR format for METIS
        xadj = (idx_t*)malloc((n + 1) * sizeof(idx_t));
        int edge_count = 0;
        xadj[0] = 0;
        
        // First pass to count edges
        for (int i = 0; i < n; i++) {
            AdjListNode* temp = graph->array[i].head;
            while (temp) {
                edge_count++;
                temp = temp->next;
            }
            xadj[i + 1] = edge_count;
        }
        
        // Allocate adjacency array
        adjncy = (idx_t*)malloc(edge_count * sizeof(idx_t));
        
        // Second pass to fill adjacency array
        edge_count = 0;
        for (int i = 0; i < n; i++) {
            AdjListNode* temp = graph->array[i].head;
            while (temp) {
                adjncy[edge_count++] = temp->dest;
                temp = temp->next;
            }
        }
        
        // Allocate partition array
        part = (idx_t*)malloc(n * sizeof(idx_t));
        
        // Call METIS to partition the graph
        int ret = METIS_PartGraphKway(&nVertices, &nConstraints, xadj, adjncy, 
                                      NULL, NULL, NULL, &nParts, NULL, 
                                      NULL, NULL, &objval, part);
                                      
        if (ret != METIS_OK) {
            printf("METIS partitioning failed with error %d\n", ret);
            exit(EXIT_FAILURE);
        }
        
        printf("Graph partitioned with edge-cut: %d\n", objval);
    }
    
    // Broadcast partition information to all processes
    int* vertex_to_proc = (int*)malloc(n * sizeof(int));
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            vertex_to_proc[i] = part[i];
        }
    }
    
    MPI_Bcast(vertex_to_proc, n, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Count vertices for this process
    int local_count = 0;
    for (int i = 0; i < n; i++) {
        if (vertex_to_proc[i] == rank) {
            local_count++;
        }
    }
    
    // Create mapping arrays
    partition->global_to_local = (int*)malloc(n * sizeof(int));
    partition->local_to_global = (int*)malloc(local_count * sizeof(int));
    
    // Set default values
    for (int i = 0; i < n; i++) {
        partition->global_to_local[i] = -1;
    }
    
    // Fill mapping arrays
    int local_idx = 0;
    for (int i = 0; i < n; i++) {
        if (vertex_to_proc[i] == rank) {
            partition->global_to_local[i] = local_idx;
            partition->local_to_global[local_idx] = i;
            local_idx++;
        }
    }
    
    // Identify boundary vertices
    partition->boundary_vertices = (int*)malloc(local_count * sizeof(int));
    partition->num_boundary = 0;
    
    for (int local_i = 0; local_i < local_count; local_i++) {
        int global_i = partition->local_to_global[local_i];
        AdjListNode* temp = graph->array[global_i].head;
        bool is_boundary = false;
        
        while (temp && !is_boundary) {
            if (vertex_to_proc[temp->dest] != rank) {
                is_boundary = true;
            }
            temp = temp->next;
        }
        
        if (is_boundary) {
            partition->boundary_vertices[partition->num_boundary++] = local_i;
        }
    }
    
    // Set partition range
    partition->start_vertex = 0;
    partition->end_vertex = local_count - 1;
    
    // Clean up
    if (rank == 0) {
        free(xadj);
        free(adjncy);
        free(part);
    }
    free(vertex_to_proc);
    
    printf("Process %d has %d vertices, %d boundary vertices\n", 
           rank, local_count, partition->num_boundary);
    
    return partition;
}

void freePartition(Partition* partition) {
    if (partition != NULL) {
        free(partition->global_to_local);
        free(partition->local_to_global);
        free(partition->boundary_vertices);
        free(partition);
    }
}

// Parallelized SOSP_Update algorithm
void SOSP_Update_Parallel(Graph* graph, SOSPTree* tree, Edge* inserted_edges, int num_inserted, 
                          int rank, int size, Partition* partition) {
    int n = graph->num_vertices;
    
    if (rank == 0) {
        printf("SOSP_Update: Parallel implementation with %d processes\n", size);
    }
    
    // Reset marked array
    for (int i = 0; i < n; i++) {
        tree->marked[i] = false;
    }
    
    // Step 0: Add edges to graph (all processes add all edges)
    for (int i = 0; i < num_inserted; i++) {
        Edge e = inserted_edges[i];
        
        // Skip invalid edges
        if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) {
            continue;
        }
        
        addEdge(graph, e.u, e.v, e.weight);
    }
    
    // Allocate arrays for affected vertices
    bool* in_affected = (bool*)calloc(n, sizeof(bool));
    int* affected = (int*)malloc(n * sizeof(int));
    int affected_size = 0;
    
    // Step 1: Process inserted edges - identify affected vertices 
    for (int i = 0; i < num_inserted; i++) {
        Edge e = inserted_edges[i];
        
        // Skip invalid edges
        if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) {
            continue;
        }
        
        int u = e.u;
        int v = e.v;
        int weight = e.weight;
        
        // Check if the distance to v can be improved
        if (tree->distance[u] != INT_MAX && 
            (tree->distance[v] == INT_MAX || tree->distance[v] > tree->distance[u] + weight)) {
            
            tree->distance[v] = tree->distance[u] + weight;
            tree->parent[v] = u;
            tree->marked[v] = true;
            
            if (!in_affected[v]) {
                in_affected[v] = true;
                affected[affected_size++] = v;
            }
        }
    }
    
    if (rank == 0) {
        printf("Initial affected vertices: %d\n", affected_size);
    }
    
    // Step 2: Propagate updates in parallel
    int* next_frontier = (int*)malloc(n * sizeof(int));
    int next_size = 0;
    bool* visited = (bool*)calloc(n, sizeof(bool));
    
    // Mark initial affected vertices as visited
    for (int i = 0; i < affected_size; i++) {
        visited[affected[i]] = true;
    }
    
    // Arrays for MPI communication
    int* all_distances = (int*)malloc(n * sizeof(int));
    int* all_parents = (int*)malloc(n * sizeof(int));
    bool* all_marked = (bool*)malloc(n * sizeof(bool));
    bool global_continue = true;
    
    // Propagation loop
    while (global_continue) {
        next_size = 0;
        
        // Process current frontier (each process handles its partition)
        for (int i = 0; i < affected_size; i++) {
            int u = affected[i];
            
            // Skip unaffected vertices
            if (!tree->marked[u]) continue;
            
            // Process each outgoing edge from u
            AdjListNode* current = graph->array[u].head;
            while (current) {
                int v = current->dest;
                int weight = current->weight;
                
                // Try to relax edge (u,v)
                if (tree->distance[u] != INT_MAX && 
                    (tree->distance[v] == INT_MAX || tree->distance[v] > tree->distance[u] + weight)) {
                    
                    tree->distance[v] = tree->distance[u] + weight;
                    tree->parent[v] = u;
                    tree->marked[v] = true;
                    
                    // Add to next frontier if not already visited
                    if (!visited[v]) {
                        visited[v] = true;
                        next_frontier[next_size++] = v;
                    }
                }
                
                current = current->next;
            }
        }
        
        // Synchronize distances, parents, and marked arrays across all processes
        MPI_Allreduce(tree->distance, all_distances, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(tree->parent, all_parents, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(tree->marked, all_marked, n, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        
        // Update local arrays with global information
        for (int i = 0; i < n; i++) {
            tree->distance[i] = all_distances[i];
            tree->parent[i] = all_parents[i];
            tree->marked[i] = all_marked[i];
        }
        
        // Determine if any process has vertices to process in the next iteration
        int local_continue = (next_size > 0);
        MPI_Allreduce(&local_continue, &global_continue, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        
        // Share next_frontier information
        int* all_next_sizes = NULL;
        int* displs = NULL;
        int* all_next_frontiers = NULL;
        
        if (rank == 0) {
            all_next_sizes = (int*)malloc(size * sizeof(int));
            displs = (int*)malloc(size * sizeof(int));
        }
        
        // Gather frontier sizes
        MPI_Gather(&next_size, 1, MPI_INT, all_next_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // Process 0 calculates displacements and total size
        int total_next_size = 0;
        if (rank == 0) {
            displs[0] = 0;
            for (int i = 0; i < size; i++) {
                total_next_size += all_next_sizes[i];
                if (i > 0) {
                    displs[i] = displs[i-1] + all_next_sizes[i-1];
                }
            }
            all_next_frontiers = (int*)malloc(total_next_size * sizeof(int));
        }
        
        // Gather all frontier vertices
        MPI_Gatherv(next_frontier, next_size, MPI_INT, all_next_frontiers, 
                   all_next_sizes, displs, MPI_INT, 0, MPI_COMM_WORLD);
        
        // Process 0 merges and broadcasts the new frontier
        if (rank == 0) {
            // Reset visited status for the new frontier
            for (int i = 0; i < n; i++) {
                visited[i] = false;
            }
            
            // Deduplicate the frontier
            affected_size = 0;
            for (int i = 0; i < total_next_size; i++) {
                int v = all_next_frontiers[i];
                if (!visited[v]) {
                    visited[v] = true;
                    affected[affected_size++] = v;
                }
            }
            
            // Reset visited for the next propagation loop
            for (int i = 0; i < n; i++) {
                visited[i] = false;
            }
            
            // Mark the new frontier as visited
            for (int i = 0; i < affected_size; i++) {
                visited[affected[i]] = true;
            }
            
            free(all_next_sizes);
            free(displs);
            free(all_next_frontiers);
        }
        
        // Broadcast the new frontier size and vertices
        MPI_Bcast(&affected_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(affected, affected_size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(visited, n, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        
        if (rank == 0 && affected_size > 0) {
            printf("Next iteration affected vertices: %d\n", affected_size);
        }
    }
    
    // Clean up
    free(affected);
    free(in_affected);
    free(next_frontier);
    free(visited);
    free(all_distances);
    free(all_parents);
    free(all_marked);
}

// Sequential Bellman-Ford implementation with parallelization
void computeInitialSOSP_Parallel(Graph* graph, SOSPTree* tree, int source, 
                                int rank, int size, Partition* partition) {
    int n = graph->num_vertices;
    
    // Initialize distances
    for (int i = 0; i < n; i++) {
        tree->distance[i] = (i == source) ? 0 : INT_MAX;
        tree->parent[i] = -1;
    }
    
    // Arrays for communication
    int* all_distances = (int*)malloc(n * sizeof(int));
    int* all_parents = (int*)malloc(n * sizeof(int));
    
    // Use a flag to check if any relaxation occurs in an iteration
    bool global_any_update = true;
    bool local_any_update = false;
    
    // Relax edges |V| - 1 times (or until no more updates)
    for (int iter = 0; iter < n - 1 && global_any_update; iter++) {
        local_any_update = false;
        
        // Each process handles its partition
        for (int idx = partition->start_vertex; idx <= partition->end_vertex; idx++) {
            int u = partition->local_to_global[idx];
            
            if (tree->distance[u] == INT_MAX) continue;
            
            AdjListNode* current = graph->array[u].head;
            while (current) {
                int v = current->dest;
                int weight = current->weight;
                
                if (tree->distance[u] != INT_MAX && 
                    (tree->distance[v] == INT_MAX || tree->distance[u] + weight < tree->distance[v])) {
                    tree->distance[v] = tree->distance[u] + weight;
                    tree->parent[v] = u;
                    local_any_update = true;
                }
                
                current = current->next;
            }
        }
        
        // Synchronize distances and parents
        MPI_Allreduce(tree->distance, all_distances, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(tree->parent, all_parents, n, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        
        // Update local copies
        memcpy(tree->distance, all_distances, n * sizeof(int));
        memcpy(tree->parent, all_parents, n * sizeof(int));
        
        // Check if any process made an update
        MPI_Allreduce(&local_any_update, &global_any_update, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        
        // Report progress for large graphs
        if (rank == 0 && (iter + 1) % 100 == 0) {
            printf("Completed %d Bellman-Ford iterations\n", iter + 1);
        }
    }
    
    free(all_distances);
    free(all_parents);
}

// Parallel MOSP_Update for 2 objectives
void MOSP_Update_2obj_Parallel(Graph* graph1, Graph* graph2, SOSPTree* tree1, SOSPTree* tree2, 
                               Edge* inserted_edges1, Edge* inserted_edges2, 
                               int num_inserted1, int num_inserted2, int source,
                               int rank, int size, Partition* partition) {
    
    int n = graph1->num_vertices;
    
    if (rank == 0) {
        printf("Step 1: Updating SOSP trees in parallel...\n");
    }
    
    // Step 1: Update SOSP trees in parallel
    double start_time_sosp1 = MPI_Wtime();
    SOSP_Update_Parallel(graph1, tree1, inserted_edges1, num_inserted1, rank, size, partition);
    double end_time_sosp1 = MPI_Wtime();
    
    double start_time_sosp2 = MPI_Wtime();
    SOSP_Update_Parallel(graph2, tree2, inserted_edges2, num_inserted2, rank, size, partition);
    double end_time_sosp2 = MPI_Wtime();
    
    if (rank == 0) {
        printf("SOSP1 update time: %.6f seconds\n", end_time_sosp1 - start_time_sosp1);
        printf("SOSP2 update time: %.6f seconds\n", end_time_sosp2 - start_time_sosp2);
        printf("Step 2: Creating combined graph...\n");
    }
    
    // Step 2: Create combined graph (each process creates its part)
    double start_time_combined = MPI_Wtime();
    Graph* combined = createGraph(n);
    
    // Add edges from both SOSP trees to combined graph
    for (int idx = partition->start_vertex; idx <= partition->end_vertex; idx++) {
        int v = partition->local_to_global[idx];
        int parent1 = tree1->parent[v];
        int parent2 = tree2->parent[v];
        
        // Add edge from parent1 to v if it exists
        if (parent1 != -1) {
            // If this edge appears in both trees, weight is 1, otherwise 2
            int weight = (parent1 == parent2) ? 1 : 2;
            addEdge(combined, parent1, v, weight);
        }
        
        // If parent2 is different from parent1, add that edge too
        if (parent2 != -1 && parent2 != parent1) {
            addEdge(combined, parent2, v, 2);
        }
    }
    
    // Barrier to ensure all processes have built their part of the combined graph
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time_combined = MPI_Wtime();
    
    if (rank == 0) {
        printf("Combined graph creation time: %.6f seconds\n", end_time_combined - start_time_combined);
        printf("Step 3: Finding SOSP in combined graph...\n");
    }
    
    // Step 3: Find SOSP in combined graph
    double start_time_final = MPI_Wtime();
    SOSPTree* combined_tree = initSOSPTree(n, source);
    
    // Use parallel Bellman-Ford to compute SOSP in combined graph
    computeInitialSOSP_Parallel(combined, combined_tree, source, rank, size, partition);
    
    double end_time_final = MPI_Wtime();
    
    if (rank == 0) {
        printf("Final SOSP computation time: %.6f seconds\n", end_time_final - start_time_final);
        
        // Find reachable vertices
        int reachable_count = 0;
        for (int i = 0; i < n; i++) {
            if (combined_tree->distance[i] != INT_MAX) {
                reachable_count++;
            }
        }
        
        printf("Vertices reachable from source in MOSP: %d\n", reachable_count);
        
        // Select a sample reachable vertex for demonstration
        int test_vertex = -1;
        for (int i = 100; i < n && test_vertex == -1; i++) {
            if (combined_tree->distance[i] != INT_MAX) {
                test_vertex = i;
                break;
            }
        }
        
        if (test_vertex != -1) {
            printf("\nMOSP result (path from source %d to vertex %d):\n", source, test_vertex);
            
            int current = test_vertex;
            int path[1000];  // Assuming path won't be longer than 1000 nodes
            int path_length = 0;
            
            // Build path in reverse
            while (current != -1 && current != source) {
                path[path_length++] = current;
                current = combined_tree->parent[current];
                
                // Safety check to avoid infinite loops
                if (path_length >= 1000) {
                    printf("Path too long or cycle detected!\n");
                    break;
                }
            }
            
            if (current == source) {
                path[path_length++] = source;
                
                // Print path in reverse order
                printf("Path: ");
                for (int i = path_length - 1; i >= 0; i--) {
                    printf("%d", path[i]);
                    if (i > 0) printf(" -> ");
                }
                printf("\n");
                
                // Print objective values along the path
                printf("Objective 1 distance: %d\n", tree1->distance[test_vertex]);
                printf("Objective 2 distance: %d\n", tree2->distance[test_vertex]);
            } else {
                printf("Error: Failed to construct path - cycle detected\n");
            }
        } else {
            printf("\nNo vertices are reachable from source vertex %d!\n", source);
        }
    }
    
    // Free memory
    freeSOSPTree(combined_tree);
    freeGraph(combined);
}

// Function to read graph from file with improved buffering
Graph* readGraph(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL;
    }
    
    char line[256];
    int num_nodes = 0, num_edges = 0;
    
    // Skip comments and find number of nodes and edges
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#') {
            // Try to extract nodes and edges info
            if (strstr(line, "Nodes:") && strstr(line, "Edges:")) {
                sscanf(line, "# Nodes: %d Edges: %d", &num_nodes, &num_edges);
                printf("Found in header: Nodes: %d, Edges: %d\n", num_nodes, num_edges);
            }
            continue;
        }
        break;
    }
    
    // If we couldn't find the node count in the header, determine max node ID
    if (num_nodes == 0) {
        rewind(file);
        
        // Skip comments
        while (fgets(line, sizeof(line), file)) {
            if (line[0] != '#') break;
        }
        
        // Find maximum node ID
        int max_node = -1;
        do {
            int u, v;
            if (sscanf(line, "%d %d", &u, &v) == 2) {
                if (u > max_node) max_node = u;
                if (v > max_node) max_node = v;
            }
        } while (fgets(line, sizeof(line), file));
        
        num_nodes = max_node + 1;  // Nodes are 0-indexed
        printf("Determined from data: Max node ID: %d\n", max_node);
    }
    
    // Create graph
    Graph* graph = createGraph(num_nodes);
    
    // Rewind and read edges
    rewind(file);
    
    // Skip comments
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '#') break;
    }
    
    // Process edges in chunks for better memory locality
    const int CHUNK_SIZE = 100000;
    Edge* edge_chunk = (Edge*)malloc(CHUNK_SIZE * sizeof(Edge));
    int chunk_count = 0;
    
    // Process the first non-comment line and continue reading
    int edges_read = 0;
    do {
        int u, v;
        if (sscanf(line, "%d %d", &u, &v) == 2) {
            // For testing purposes, assign random weights
            int weight = rand() % 10 + 1;
            
            // Check if nodes are within bounds
            if (u >= 0 && u < num_nodes && v >= 0 && v < num_nodes) {
                edge_chunk[chunk_count].u = u;
                edge_chunk[chunk_count].v = v;
                edge_chunk[chunk_count].weight = weight;
                chunk_count++;
                edges_read++;
                
                // Process chunk when it's full
                if (chunk_count == CHUNK_SIZE) {
                    for (int i = 0; i < chunk_count; i++) {
                        addEdge(graph, edge_chunk[i].u, edge_chunk[i].v, edge_chunk[i].weight);
                    }
                    chunk_count = 0;
                }
                
                // Progress indicator for large files
                if (edges_read % 1000000 == 0) {
                    printf("Read %d edges...\n", edges_read);
                }
            }
        }
    } while (fgets(line, sizeof(line), file));
    
    // Process remaining edges in the last chunk
    for (int i = 0; i < chunk_count; i++) {
        addEdge(graph, edge_chunk[i].u, edge_chunk[i].v, edge_chunk[i].weight);
    }
    
    free(edge_chunk);
    fclose(file);
    printf("Graph created with %d vertices and %d edges\n", graph->num_vertices, graph->num_edges);
    
    return graph;
}

// Add this sequential version if it's not already defined
void computeInitialSOSP(Graph* graph, SOSPTree* tree, int source) {
    int n = graph->num_vertices;
    
    // Initialize distances
    for (int i = 0; i < n; i++) {
        tree->distance[i] = (i == source) ? 0 : INT_MAX;
        tree->parent[i] = -1;
    }
    
    // Use a flag to check if any relaxation occurs in an iteration
    bool any_update;
    
    // Relax edges |V| - 1 times (or until no more updates)
    for (int iter = 0; iter < n - 1; iter++) {
        any_update = false;
        
        for (int u = 0; u < n; u++) {
            if (tree->distance[u] == INT_MAX) continue;
            
            AdjListNode* current = graph->array[u].head;
            while (current) {
                int v = current->dest;
                int weight = current->weight;
                
                if (tree->distance[u] + weight < tree->distance[v]) {
                    tree->distance[v] = tree->distance[u] + weight;
                    tree->parent[v] = u;
                    any_update = true;
                }
                
                current = current->next;
            }
        }
        
        // Stop if no more relaxations
        if (!any_update) {
            printf("Bellman-Ford converged after %d iterations\n", iter + 1);
            break;
        }
        
        // Report progress for large graphs
        if ((iter + 1) % 100 == 0) {
            printf("Completed %d Bellman-Ford iterations\n", iter + 1);
        }
    }
}

// Modified createSmartRandomEdges function using sequential computeInitialSOSP
Edge* createSmartRandomEdges(Graph* graph, int source, int target, int num_edges) {
    int n = graph->num_vertices;
    Edge* edges = (Edge*)malloc(num_edges * sizeof(Edge));
    if (edges == NULL) {
        fprintf(stderr, "Memory allocation failed for random edges\n");
        exit(EXIT_FAILURE);
    }
    
    // Ensure some edges connect the source to target components
    int connected_edges = num_edges / 4;  // 25% of edges will help connectivity
    
    // Create a temporary SOSP tree to find reachable vertices from source
    SOSPTree* temp_tree = initSOSPTree(n, source);
    computeInitialSOSP(graph, temp_tree, source);  // Use sequential version
    
    // Find vertices reachable from source
    int* reachable = (int*)malloc(n * sizeof(int));
    int reachable_count = 0;
    
    for (int i = 0; i < n; i++) {
        if (temp_tree->distance[i] != INT_MAX) {
            reachable[reachable_count++] = i;
        }
    }
    
    // Check if target is reachable
    bool target_reachable = false;
    for (int i = 0; i < reachable_count; i++) {
        if (reachable[i] == target) {
            target_reachable = true;
            break;
        }
    }
    
    printf("Initial analysis: %d vertices reachable from source. Target %s reachable.\n", 
           reachable_count, target_reachable ? "is" : "is NOT");
    
    // If target is not reachable, add edges to make it reachable
    if (!target_reachable && reachable_count > 0) {
        // Find vertices that can reach target
        freeSOSPTree(temp_tree);
        temp_tree = initSOSPTree(n, target);
        
        // Create a reversed graph to find vertices that can reach target
        Graph* reversed = createGraph(n);
        for (int u = 0; u < n; u++) {
            AdjListNode* current = graph->array[u].head;
            while (current) {
                int v = current->dest;
                addEdge(reversed, v, u, current->weight);
                current = current->next;
            }
        }
        
        computeInitialSOSP(reversed, temp_tree, target);  // Use sequential version
        
        int* can_reach_target = (int*)malloc(n * sizeof(int));
        int can_reach_count = 0;
        for (int i = 0; i < n; i++) {
            if (temp_tree->distance[i] != INT_MAX) {
                can_reach_target[can_reach_count++] = i;
            }
        }
        
        printf("Found %d vertices that can reach target.\n", can_reach_count);
        
        // Create edges that connect reachable_from_source to can_reach_target
        for (int i = 0; i < connected_edges && i < reachable_count && can_reach_count > 0; i++) {
            int u = reachable[rand() % reachable_count];
            int v = can_reach_target[rand() % can_reach_count];
            
            edges[i].u = u;
            edges[i].v = v;
            edges[i].weight = rand() % 10 + 1;
            
            printf("Creating connecting edge: %d -> %d\n", u, v);
        }
        
        free(can_reach_target);
        freeGraph(reversed);
    } else {
        // If target is already reachable, create some connecting edges anyway
        for (int i = 0; i < connected_edges && reachable_count > 0; i++) {
            edges[i].u = reachable[rand() % reachable_count];
            edges[i].v = rand() % n;
            edges[i].weight = rand() % 10 + 1;
        }
    }
    
    // Fill the rest with random edges
    for (int i = connected_edges; i < num_edges; i++) {
        if (reachable_count > 0) {
            // 50% of edges start from reachable vertices
            if (rand() % 2 == 0) {
                edges[i].u = reachable[rand() % reachable_count];
            } else {
                edges[i].u = rand() % n;
            }
        } else {
            edges[i].u = rand() % n;
        }
        
        edges[i].v = rand() % n;
        edges[i].weight = rand() % 10 + 1;
    }
    
    free(reachable);
    freeSOSPTree(temp_tree);
    
    return edges;
}

// Main function with MPI initialization
int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: %s <graph_file> [source_vertex] [target_vertex]\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    // Record total execution time
    double total_start_time = MPI_Wtime();
    
    // Set random seed based on rank to avoid identical random numbers
    srand(time(NULL) + rank);
    
    const char* filename = argv[1];
    int source_vertex = 0;  // Default source vertex
    int target_vertex = 100; // Default target vertex
    
    // Override defaults if provided
    if (argc >= 3) {
        source_vertex = atoi(argv[2]);
    }
    
    if (argc >= 4) {
        target_vertex = atoi(argv[3]);
    }
    
    if (rank == 0) {
        printf("Parallel MOSP Update Implementation with MPI (%d processes)\n", size);
    }
    
    // Read graphs for two objectives (all processes read the graph)
    if (rank == 0) {
        printf("Reading graph from file %s...\n", filename);
    }
    double start_time_read = MPI_Wtime();
    Graph* graph1 = readGraph(filename);
    Graph* graph2 = readGraph(filename);  // Same graph but different weights will be used
    double end_time_read = MPI_Wtime();
    
    if (rank == 0) {
        printf("Graph reading time: %.6f seconds\n", end_time_read - start_time_read);
    }
    
    if (!graph1 || !graph2) {
        if (rank == 0) {
            printf("Failed to read graph from file\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    int n = graph1->num_vertices;
    
    if (source_vertex >= n || source_vertex < 0) {
        if (rank == 0) {
            printf("Source vertex %d is out of range (0-%d). Using default source 0.\n", 
                   source_vertex, n-1);
        }
        source_vertex = 0;
    }
    
    if (target_vertex >= n || target_vertex < 0) {
        if (rank == 0) {
            printf("Target vertex %d is out of range (0-%d). Using vertex %d instead.\n", 
                   target_vertex, n-1, n/2);
        }
        target_vertex = n/2;
    }
    
    // Partition the graph using METIS
    if (rank == 0) {
        printf("Partitioning graph for %d processes...\n", size);
    }
    double start_time_partition = MPI_Wtime();
    Partition* partition = partitionGraphMETIS(graph1, rank, size);
    double end_time_partition = MPI_Wtime();
    
    if (rank == 0) {
        printf("Graph partitioning time: %.6f seconds\n", 
               end_time_partition - start_time_partition);
    }
    
    // Initialize SOSP trees
    if (rank == 0) {
        printf("Initializing SOSP trees for source vertex %d...\n", source_vertex);
    }
    SOSPTree* tree1 = initSOSPTree(n, source_vertex);
    SOSPTree* tree2 = initSOSPTree(n, source_vertex);
    
    // Compute initial SOSP trees in parallel
    if (rank == 0) {
        printf("Computing initial SOSP trees in parallel...\n");
    }
    double start_time_initial = MPI_Wtime();
    computeInitialSOSP_Parallel(graph1, tree1, source_vertex, rank, size, partition);
    computeInitialSOSP_Parallel(graph2, tree2, source_vertex, rank, size, partition);
    double end_time_initial = MPI_Wtime();
    
    if (rank == 0) {
        printf("Initial SOSP computation time: %.6f seconds\n", 
               end_time_initial - start_time_initial);
    }
    
    // Generate inserted edges (only process 0 generates them)
    int num_inserted = 5000;  // Number of edges to insert
    Edge* inserted_edges1 = NULL;
    Edge* inserted_edges2 = NULL;
    
    if (rank == 0) {
        printf("Generating %d 'smart' random edges to insert...\n", num_inserted);
        
        double start_time_edges = MPI_Wtime();
        inserted_edges1 = createSmartRandomEdges(graph1, source_vertex, target_vertex, num_inserted);
        inserted_edges2 = createSmartRandomEdges(graph2, source_vertex, target_vertex, num_inserted);
        double end_time_edges = MPI_Wtime();
        
        printf("Edge generation time: %.6f seconds\n", 
               end_time_edges - start_time_edges);
    }
    
    // Broadcast inserted edges to all processes
    if (inserted_edges1 == NULL) {
        inserted_edges1 = (Edge*)malloc(num_inserted * sizeof(Edge));
    }
    if (inserted_edges2 == NULL) {
        inserted_edges2 = (Edge*)malloc(num_inserted * sizeof(Edge));
    }
    
    MPI_Bcast(inserted_edges1, num_inserted * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(inserted_edges2, num_inserted * sizeof(Edge), MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Run MOSP_Update in parallel
    if (rank == 0) {
        printf("Running MOSP_Update in parallel...\n");
    }
    double start_time = MPI_Wtime();
    
    MOSP_Update_2obj_Parallel(graph1, graph2, tree1, tree2, 
                             inserted_edges1, inserted_edges2, 
                             num_inserted, num_inserted, source_vertex,
                             rank, size, partition);
    
    double end_time = MPI_Wtime();
    
    if (rank == 0) {
        printf("\nMOSP_Update completed in %.6f seconds\n", end_time - start_time);
    }
    
    // Clean up
    freeSOSPTree(tree1);
    freeSOSPTree(tree2);
    freeGraph(graph1);
    freeGraph(graph2);
    free(inserted_edges1);
    free(inserted_edges2);
    freePartition(partition);
    
    // Calculate and print total execution time
    double total_end_time = MPI_Wtime();
    
    if (rank == 0) {
        printf("\nTotal program execution time: %.6f seconds\n", 
               total_end_time - total_start_time);
    }
    
    MPI_Finalize();
    return 0;
}
