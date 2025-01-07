#include "binary_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <time.h>
#include <float.h>
#include <math.h>

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define CYAN    "\x1b[36m"    
#define RESET   "\x1b[0m"

#define TEST_ITERATIONS 10
#define NUM_THREADS 8
#define STRESS_SIZE 1000

// Add statistics structure
typedef struct {
    double min_time;
    double max_time;
    double total_time;
    int successful_runs;
} TestStats;

// Add statistics initialization
static void init_stats(TestStats* stats) {
    stats->min_time = DBL_MAX;
    stats->max_time = 0;
    stats->total_time = 0;
    stats->successful_runs = 0;
}

// Add statistics reporting
static void report_stats(const char* test_name, TestStats* stats) {
    printf(CYAN "\nStatistics for %s (%d successful runs):\n" RESET, 
           test_name, stats->successful_runs);
    printf(YELLOW "Average time: %f seconds\n" RESET, 
           stats->total_time / stats->successful_runs);
    printf(YELLOW "Min time: %f seconds\n" RESET, stats->min_time);
    printf(YELLOW "Max time: %f seconds\n" RESET, stats->max_time);
}

// Independent verification functions
static bool verify_search(TreeNode* root, int data, bool expected) {
    if (root == NULL) {
        return !expected;
    }
    
    TreeNode* current = root;
    while (current != NULL) {
        if (current->data == data) {
            return expected;
        }
        current = data < current->data ? current->left : current->right;
    }
    return !expected;
}

static bool verify_tree_empty(TreeNode* root) {
    return root == NULL;
}

static bool verify_node_exists(TreeNode* root, int data) {
    if (root == NULL) {
        return false;
    }
    
    TreeNode* current = root;
    while (current != NULL) {
        if (current->data == data) {
            return true;
        }
        if (data < current->data) {
            current = current->left;
        } else {
            current = current->right;
        }
    }
    return false;
}

static bool verify_binary_search_property(TreeNode* root) {
    if (root == NULL) {
        return true;
    }
    
    // Check if left subtree contains only smaller values
    TreeNode* left = root->left;
    while (left != NULL) {
        if (left->data > root->data) {
            printf(RED "VIOLATION: %d >= %d\n" RESET, left->data, root->data);
            return false;
        }
        left = left->right;  // Check rightmost node of left subtree
    }
    
    // Check if right subtree contains only larger values
    TreeNode* right = root->right;
    while (right != NULL) {
        if (right->data < root->data) {
            printf(RED "VIOLATION: %d <= %d\n" RESET, right->data, root->data);
            return false;
        }
        right = right->left;  // Check leftmost node of right subtree
    }
    
    // Recursively check subtrees
    return verify_binary_search_property(root->left) && 
           verify_binary_search_property(root->right);
}

static void print_tree_structure(TreeNode* root, int level) {
    if (root == NULL) {
        return;
    }
    
    print_tree_structure(root->right, level + 1);
    
    // Print indentation and node
    printf("\033[0m");  // Reset color
    for(int i = 0; i < level; i++) {
        printf("    ");  // 4 spaces for each level
    }
    
    // Color-code based on level
    if (level == 0) {
        printf(CYAN "┌─[%d]" RESET "\n", root->data);  // Root node
    } else {
        if (level % 3 == 0) {
            printf(YELLOW "├─[%d]" RESET "\n", root->data);
        } else if (level % 3 == 1) {
            printf(GREEN "├─[%d]" RESET "\n", root->data);
        } else {
            printf(BLUE "├─[%d]" RESET "\n", root->data);
        }
    }
    
    print_tree_structure(root->left, level + 1);
}

// Add this to test functions to verify tree structure
static void verify_tree_structure(TreeNode* root) {
    if (!verify_binary_search_property(root)) {
        printf(RED "ERROR: Binary search tree property violated\n" RESET);
        printf("Tree structure:\n");
        print_tree_structure(root, 0);
        exit(1);
    }
}

void test_basic_operations() {
    printf(BLUE "Testing basic operations...\n" RESET);
    TreeNode* root = NULL;
    
    // Test complex insertion pattern
    int values[] = {50, 30, 70, 20, 40, 60, 80, 15, 25, 35, 45, 55, 65, 75, 85};
    int n = sizeof(values) / sizeof(values[0]);
    
    for(int i = 0; i < n; i++) {
        root = insertNode(root, values[i]);
        verify_tree_structure(root);
    }
    
    // Test search for all inserted values
    for(int i = 0; i < n; i++) {
        if (!verify_node_exists(root, values[i])) {
            printf(RED "ERROR: Node %d not found after insertion\n" RESET, values[i]);
            exit(1);
        }
    }
    
    // Test deletion with different node types (leaf, one child, two children)
    root = deleteNode(root, 15);  // leaf node
    root = deleteNode(root, 60);  // one child
    root = deleteNode(root, 30);  // two children
    
    if (!verify_search(root, 15, false) || 
        !verify_search(root, 60, false) || 
        !verify_search(root, 30, false)) {
        printf(RED "ERROR: Deletion test failed\n" RESET);
        exit(1);
    }
    
    freeTree(root);
    printf(GREEN "Basic operations test passed!\n" RESET);
}

void test_parallel_operations() {
    printf(BLUE "Testing parallel operations...\n" RESET);
    TreeNode* root = NULL;
    root = insertNode(root, STRESS_SIZE/2);  // Start with middle value
    
    double start_time = omp_get_wtime();
    
    // Parallel insertion
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for nowait
        for (int i = 0; i < STRESS_SIZE; i++) {
            if (i != STRESS_SIZE/2) { // Skip the root value
                insertNode(root, i);
            }
        }
    }
    
    double insert_time = omp_get_wtime() - start_time;
    printf(YELLOW "Parallel insertion time: %f seconds\n" RESET, insert_time);
    
    // Verify all insertions
    bool all_found = true;
    #pragma omp parallel num_threads(NUM_THREADS) reduction(&:all_found)
    {
        #pragma omp for nowait
        for (int i = 0; i < STRESS_SIZE; i++) {
            bool exists = verify_node_exists(root, i);
            all_found &= exists;
            if (!exists) {
                printf(RED "ERROR: Node %d not found after insertion\n" RESET, i);
            }
        }
    }
    
    if (!all_found) {
        printf(RED "ERROR: Not all nodes were inserted correctly\n" RESET);
        exit(1);
    }
    
    // Parallel search
    start_time = omp_get_wtime();
    int found_count = 0;
    
    #pragma omp parallel num_threads(NUM_THREADS) reduction(+:found_count)
    {
        #pragma omp for nowait
        for (int i = 0; i < STRESS_SIZE; i++) {
            if (searchNode(root, i)) {
                found_count++;
            }
        }
    }
    
    double search_time = omp_get_wtime() - start_time;
    printf(YELLOW "Parallel search time: %f seconds\n" RESET, search_time);
    assert(found_count == STRESS_SIZE);
    
    // Parallel deletion
    start_time = omp_get_wtime();
    
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for nowait
        for (int i = 0; i < STRESS_SIZE/2; i++) {
            deleteNode(root, i);
        }
    }
    
    double delete_time = omp_get_wtime() - start_time;
    printf(YELLOW "Parallel deletion time: %f seconds\n" RESET, delete_time);
    
    // Add structure verification after parallel operations
    verify_tree_structure(root);
    
    freeTree(root);
    printf(GREEN "Parallel operations test passed!\n" RESET);
}

void test_traversals() {
    printf("Testing traversals...\n");
    TreeNode* root = NULL;
    
    root = insertNode(root, 50);
    root = insertNode(root, 30);
    root = insertNode(root, 70);
    root = insertNode(root, 20);
    root = insertNode(root, 40);
    
    printf("Inorder traversal: ");
    inorderTraversal(root);
    printf("\n");
    
    printf("Preorder traversal: ");
    preorderTraversal(root);
    printf("\n");
    
    printf("Postorder traversal: ");
    postorderTraversal(root);
    printf("\n");
    
    freeTree(root);
    printf("Traversal tests completed!\n");
}

void test_edge_cases() {
    printf(BLUE "Testing edge cases...\n" RESET);
    TreeNode* root = NULL;
    
    if (!verify_search(root, 10, false)) {
        printf(RED "ERROR: Empty tree search test failed\n" RESET);
        exit(1);
    }
    
    root = deleteNode(root, 10);
    if (!verify_tree_empty(root)) {
        printf(RED "ERROR: Empty tree deletion test failed\n" RESET);
        exit(1);
    }
    
    root = insertNode(root, 10);
    if (!verify_node_exists(root, 10)) {
        printf(RED "ERROR: Single node insertion test failed\n" RESET);
        exit(1);
    }
    
    root = deleteNode(root, 10);
    if (!verify_tree_empty(root)) {
        printf(RED "ERROR: Single node deletion test failed\n" RESET);
        exit(1);
    }
    
    printf(GREEN "Edge cases test passed!\n" RESET);
}

void test_unbalanced_tree() {
    printf(BLUE "Testing unbalanced tree operations...\n" RESET);
    TreeNode* root = NULL;
    
    // Create a deliberately unbalanced tree
    for(int i = 0; i < 1000; i++) {
        root = insertNode(root, i);
    }
    
    double start_time = omp_get_wtime();
    
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for nowait
        for(int i = 0; i < 1000; i++) {
            assert(searchNode(root, i) == true);
        }
    }
    
    double search_time = omp_get_wtime() - start_time;
    printf(YELLOW "Unbalanced tree search time: %f seconds\n" RESET, search_time);
    
    freeTree(root);
    printf(GREEN "Unbalanced tree test passed!\n" RESET);
}

void test_random_operations() {
    printf(BLUE "Testing random mixed operations...\n" RESET);
    TreeNode* root = NULL;
    root = insertNode(root, 500);  // Dummy node to avoid empty tree
    srand(time(NULL));
    
    // Track operations for verification
    int inserted_values[1000] = {0};
    int insert_count = 0;
    int delete_count = 0;
    
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for nowait
        for(int i = 0; i < STRESS_SIZE; i++) {
            int op = rand() % 3;  // 0: insert, 1: search, 2: delete
            int value = rand() % 1000;
            
            switch(op) {
                case 0:
                    insertNode(root, value);
                    #pragma omp critical
                    {
                        inserted_values[insert_count++] = value;
                    }
                    break;
                case 1:
                    searchNode(root, value);
                    break;
                case 2:
                    deleteNode(root, value);
                    #pragma omp critical
                    {
                        delete_count++;
                    }
                    break;
            }
        }
    }
    
    // Verify tree structure
    verify_tree_structure(root);
    
    // Verify that all remaining inserted values are searchable
    int found_count = 0;
    for(int i = 0; i < insert_count; i++) {
        if (searchNode(root, inserted_values[i])) {
            found_count++;
        }
    }
    
    printf(YELLOW "Operations summary:\n");
    printf("- Insertions: %d\n", insert_count);
    printf("- Deletions: %d\n", delete_count);
    printf("- Found nodes after operations: %d\n" RESET, found_count);
    
    if (!verify_binary_search_property(root)) {
        printf(RED "ERROR: Binary search tree property violated after random operations\n" RESET);
        exit(1);
    }
    
    freeTree(root);
    printf(GREEN "Random operations test passed!\n" RESET);
}

// Modify main function to run multiple iterations
int main() {
    printf(YELLOW "Starting enhanced binary tree tests (%d iterations)...\n\n" RESET, 
           TEST_ITERATIONS);
    
    TestStats basic_stats, parallel_stats, unbalanced_stats, random_stats;
    init_stats(&basic_stats);
    init_stats(&parallel_stats);
    init_stats(&unbalanced_stats);
    init_stats(&random_stats);
    
    for(int i = 0; i < TEST_ITERATIONS; i++) {
        printf(CYAN "\nIteration %d/%d\n" RESET, i + 1, TEST_ITERATIONS);
        
        double start_time = omp_get_wtime();
        test_basic_operations();
        double end_time = omp_get_wtime();
        basic_stats.min_time = fmin(basic_stats.min_time, end_time - start_time);
        basic_stats.max_time = fmax(basic_stats.max_time, end_time - start_time);
        basic_stats.total_time += end_time - start_time;
        basic_stats.successful_runs++;
        
        start_time = omp_get_wtime();
        test_parallel_operations();
        end_time = omp_get_wtime();
        parallel_stats.min_time = fmin(parallel_stats.min_time, end_time - start_time);
        parallel_stats.max_time = fmax(parallel_stats.max_time, end_time - start_time);
        parallel_stats.total_time += end_time - start_time;
        parallel_stats.successful_runs++;
        
        start_time = omp_get_wtime();
        test_unbalanced_tree();
        end_time = omp_get_wtime();
        unbalanced_stats.min_time = fmin(unbalanced_stats.min_time, end_time - start_time);
        unbalanced_stats.max_time = fmax(unbalanced_stats.max_time, end_time - start_time);
        unbalanced_stats.total_time += end_time - start_time;
        unbalanced_stats.successful_runs++;
        
        start_time = omp_get_wtime();
        test_random_operations();
        end_time = omp_get_wtime();
        random_stats.min_time = fmin(random_stats.min_time, end_time - start_time);
        random_stats.max_time = fmax(random_stats.max_time, end_time - start_time);
        random_stats.total_time += end_time - start_time;
        random_stats.successful_runs++;
        
        printf(GREEN "Iteration %d completed successfully!\n" RESET, i + 1);
    }
    
    // Report final statistics
    report_stats("Basic Operations", &basic_stats);
    report_stats("Parallel Operations", &parallel_stats);
    report_stats("Unbalanced Tree", &unbalanced_stats);
    report_stats("Random Operations", &random_stats);
    
    printf(GREEN "\nAll test iterations completed successfully!\n" RESET);
    return 0;
}
