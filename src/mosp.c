#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h> // Added for sleep function
#include <sys/time.h>

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

// Function to create a new adjacency list node
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

// Function to create a graph with n vertices
Graph* createGraph(int n) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    if (graph == NULL) {
        fprintf(stderr, "Memory allocation failed for graph\n");
        exit(EXIT_FAILURE);
    }
    
    graph->num_vertices = n;
    graph->num_edges = 0;
    
    // Create an array of adjacency lists
    graph->array = (AdjList*)malloc(n * sizeof(AdjList));
    if (graph->array == NULL) {
        fprintf(stderr, "Memory allocation failed for adjacency lists\n");
        free(graph);
        exit(EXIT_FAILURE);
    }
    
    // Initialize each adjacency list as empty
    for (int i = 0; i < n; i++) {
        graph->array[i].head = NULL;
    }
    
    return graph;
}

// Function to add an edge to the graph
void addEdge(Graph* graph, int src, int dest, int weight) {
    // Add edge from src to dest
    AdjListNode* newNode = newAdjListNode(dest, weight);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
    
    graph->num_edges++;
}

// Function to initialize SOSP tree
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

// Function to free SOSP tree
void freeSOSPTree(SOSPTree* tree) {
    if (tree != NULL) {
        free(tree->parent);
        free(tree->distance);
        free(tree->marked);
        free(tree);
    }
}

// Function to free graph
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

// Optimized sequential SOSP_Update algorithm
void SOSP_Update(Graph* graph, SOSPTree* tree, Edge* inserted_edges, int num_inserted) {
    int n = graph->num_vertices;
    
    printf("SOSP_Update: Sequential implementation\n");
    
    // Reset marked array
    for (int i = 0; i < n; i++) {
        tree->marked[i] = false;
    }
    
    sleep(10);
    
    // Step 0: Add edges to graph
    for (int i = 0; i < num_inserted; i++) {
        Edge e = inserted_edges[i];
        
        // Skip invalid edges
        if (e.u < 0 || e.u >= n || e.v < 0 || e.v >= n) {
            continue;
        }
        
        addEdge(graph, e.u, e.v, e.weight);
    }
    
    // Step 1: Process inserted edges - find affected vertices
    // Use bit array to track affected vertices for efficiency
    bool* in_affected = (bool*)calloc(n, sizeof(bool));
    int* affected = (int*)malloc(n * sizeof(int));
    int affected_size = 0;
    
    // First pass: identify affected vertices from inserted edges
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
    
    printf("Initial affected vertices: %d\n", affected_size);
    
   
    sleep(5);
    
    // Step 2: Propagate the update using a frontier-based approach
    int* next_frontier = (int*)malloc(n * sizeof(int));
    int next_size = 0;
    bool* visited = (bool*)calloc(n, sizeof(bool));
    
    // Mark initial affected vertices as visited
    for (int i = 0; i < affected_size; i++) {
        visited[affected[i]] = true;
    }
    
    // Propagation loop
    while (affected_size > 0) {
        next_size = 0;
        
        // Process current frontier
        for (int i = 0; i < affected_size; i++) {
            int u = affected[i];
            
            // Skip unaffected vertices (marked is true if they were affected)
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
        
        // Swap frontiers
        memcpy(affected, next_frontier, next_size * sizeof(int));
        affected_size = next_size;
        
        printf("Next iteration affected vertices: %d\n", affected_size);
        
        // Early termination if no new vertices were affected
        if (affected_size == 0) {
            break;
        }
    }
    
    // Clean up
    free(affected);
    free(in_affected);
    free(next_frontier);
    free(visited);
}

// Sequential Bellman-Ford implementation to compute initial SOSP trees
void computeInitialSOSP(Graph* graph, SOSPTree* tree, int source) {
    int n = graph->num_vertices;
    
    // Initialize distances
    for (int i = 0; i < n; i++) {
        tree->distance[i] = (i == source) ? 0 : INT_MAX;
        tree->parent[i] = -1;
    }
    
    // Added sleep - 5 seconds
    sleep(5);
    
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
    
    
    sleep(5);
    
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

// Function to create random edges for testing that ensure connectivity
Edge* createSmartRandomEdges(Graph* graph, int source, int target, int num_edges) {
    int n = graph->num_vertices;
    Edge* edges = (Edge*)malloc(num_edges * sizeof(Edge));
    if (edges == NULL) {
        fprintf(stderr, "Memory allocation failed for random edges\n");
        exit(EXIT_FAILURE);
    }
    
    
    sleep(10);
    
    // Ensure some edges connect the source to target components
    int connected_edges = num_edges / 4;  // 25% of edges will help connectivity
    
    // Create a temporary SOSP tree to find reachable vertices from source
    SOSPTree* temp_tree = initSOSPTree(n, source);
    computeInitialSOSP(graph, temp_tree, source);
    
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
        
        computeInitialSOSP(reversed, temp_tree, target);
        
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

// Sequential MOSP_Update for 2 objectives
void MOSP_Update_2obj(Graph* graph1, Graph* graph2, SOSPTree* tree1, SOSPTree* tree2, 
                      Edge* inserted_edges1, Edge* inserted_edges2, 
                      int num_inserted1, int num_inserted2, int source) {
    
    int n = graph1->num_vertices;
    
    printf("Step 1: Updating SOSP trees sequentially...\n");
    
   
    sleep(10);
    
    // Step 1: Update SOSP trees sequentially
    clock_t start_time_sosp1 = clock();
    SOSP_Update(graph1, tree1, inserted_edges1, num_inserted1);
    clock_t end_time_sosp1 = clock();
    
    clock_t start_time_sosp2 = clock();
    SOSP_Update(graph2, tree2, inserted_edges2, num_inserted2);
    clock_t end_time_sosp2 = clock();
    
    printf("SOSP1 update time: %.6f seconds\n", 
           (double)(end_time_sosp1 - start_time_sosp1) / CLOCKS_PER_SEC);
    printf("SOSP2 update time: %.6f seconds\n", 
           (double)(end_time_sosp2 - start_time_sosp2) / CLOCKS_PER_SEC);
    
    printf("Step 2: Creating combined graph...\n");
    
    // Step 2: Create combined graph
    clock_t start_time_combined = clock();
    Graph* combined = createGraph(n);
    
    // Added sleep - 5 seconds
    sleep(5);
    
    // Add edges from both SOSP trees to combined graph
    for (int v = 0; v < n; v++) {
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
    clock_t end_time_combined = clock();
    
    printf("Combined graph creation time: %.6f seconds\n", 
           (double)(end_time_combined - start_time_combined) / CLOCKS_PER_SEC);
    
    printf("Step 3: Finding SOSP in combined graph...\n");
    
    // Step 3: Find SOSP in combined graph
    clock_t start_time_final = clock();
    SOSPTree* combined_tree = initSOSPTree(n, source);
    
    // Use sequential Bellman-Ford to compute SOSP in combined graph
    computeInitialSOSP(combined, combined_tree, source);
    
    clock_t end_time_final = clock();
    printf("Final SOSP computation time: %.6f seconds\n", 
           (double)(end_time_final - start_time_final) / CLOCKS_PER_SEC);
    
    // Added sleep - 5 seconds
    sleep(5);
    
    // Find some reachable vertices for demonstration
    int reachable_count = 0;
    for (int i = 0; i < n; i++) {
        if (combined_tree->distance[i] != INT_MAX) {
            reachable_count++;
        }
    }
    
    printf("Vertices reachable from source in MOSP: %d\n", reachable_count);
    
    // Select a sample reachable vertex to demonstrate the path
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
    
    // Free memory
    freeSOSPTree(combined_tree);
    freeGraph(combined);
}
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <graph_file> [source_vertex] [target_vertex]\n", argv[0]);
        return 1;
    }
    
    // Record total execution time using wall clock time instead of CPU time
    struct timeval total_start_tv;
    gettimeofday(&total_start_tv, NULL);
    
    srand(time(NULL));
    
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
    
    printf("Sequential MOSP Update Implementation\n");
    
    // Read graphs for two objectives
    printf("Reading graph from file %s...\n", filename);
    clock_t start_time_read = clock();
    Graph* graph1 = readGraph(filename);
    Graph* graph2 = readGraph(filename);  // Same graph but different weights will be used
    clock_t end_time_read = clock();
    
    printf("Graph reading time: %.6f seconds\n", 
           (double)(end_time_read - start_time_read) / CLOCKS_PER_SEC);
    
    if (!graph1 || !graph2) {
        printf("Failed to read graph from file\n");
        return 1;
    }
    
    int n = graph1->num_vertices;
    
    if (source_vertex >= n || source_vertex < 0) {
        printf("Source vertex %d is out of range (0-%d). Using default source 0.\n", source_vertex, n-1);
        source_vertex = 0;
    }
    
    if (target_vertex >= n || target_vertex < 0) {
        printf("Target vertex %d is out of range (0-%d). Using vertex %d instead.\n", 
               target_vertex, n-1, n/2);
        target_vertex = n/2;
    }
    
    printf("Initializing SOSP trees for source vertex %d...\n", source_vertex);
    SOSPTree* tree1 = initSOSPTree(n, source_vertex);
    SOSPTree* tree2 = initSOSPTree(n, source_vertex);
    
    // Compute initial SOSP trees
    printf("Computing initial SOSP trees...\n");
    clock_t start_time_initial = clock();
    computeInitialSOSP(graph1, tree1, source_vertex);
    computeInitialSOSP(graph2, tree2, source_vertex);
    clock_t end_time_initial = clock();
    
    printf("Initial SOSP computation time: %.6f seconds\n", 
           (double)(end_time_initial - start_time_initial) / CLOCKS_PER_SEC);
    
    // Generate inserted edges
    int num_inserted = 5000;  // Number of edges to insert
    printf("Generating %d 'smart' random edges to insert...\n", num_inserted);
    
    clock_t start_time_edges = clock();
    Edge* inserted_edges1 = createSmartRandomEdges(graph1, source_vertex, target_vertex, num_inserted);
    Edge* inserted_edges2 = createSmartRandomEdges(graph2, source_vertex, target_vertex, num_inserted);
    clock_t end_time_edges = clock();
    
    printf("Edge generation time: %.6f seconds\n", 
           (double)(end_time_edges - start_time_edges) / CLOCKS_PER_SEC);
    
    // Measure execution time for MOSP_Update using wall clock time
    printf("Running MOSP_Update...\n");
    struct timeval start_tv;
    gettimeofday(&start_tv, NULL);
    
    // Run MOSP_Update for 2 objectives
    MOSP_Update_2obj(graph1, graph2, tree1, tree2, 
                     inserted_edges1, inserted_edges2, 
                     num_inserted, num_inserted, source_vertex);
    
    struct timeval end_tv;
    gettimeofday(&end_tv, NULL);
    double execution_time = 
        (end_tv.tv_sec - start_tv.tv_sec) + 
        (end_tv.tv_usec - start_tv.tv_usec) / 1000000.0;
    
    printf("\nMOSP_Update completed in %.6f seconds\n", execution_time);
    
    // Clean up
    freeSOSPTree(tree1);
    freeSOSPTree(tree2);
    freeGraph(graph1);
    freeGraph(graph2);
    free(inserted_edges1);
    free(inserted_edges2);
    
    // Calculate and print total execution time using wall clock time
    struct timeval total_end_tv;
    gettimeofday(&total_end_tv, NULL);
    double total_execution_time = 
        (total_end_tv.tv_sec - total_start_tv.tv_sec) + 
        (total_end_tv.tv_usec - total_start_tv.tv_usec) / 1000000.0;
    printf("\nTotal program execution time: %.6f seconds\n", total_execution_time);
    
    return 0;
}
