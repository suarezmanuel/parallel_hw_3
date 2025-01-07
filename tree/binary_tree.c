#include "binary_tree.h"
#include <stdlib.h>
#include <stdio.h>

TreeNode* createNode(int data) {
    TreeNode* node = (TreeNode*) malloc(sizeof(TreeNode));
    if (node == NULL) { perror("allocation error"); return NULL; }
    node->data = data;
    omp_init_lock(&node->lock);
    node->left = NULL;
    node->right = NULL;
    return node;
}

TreeNode* insertNode(TreeNode* root, int data) {

    if (root == NULL) return createNode(data);

    TreeNode* curr = root;
    omp_set_lock(&curr->lock);

    TreeNode* next = (data <= curr->data) ? curr->left : curr->right;

    while (next != NULL) {

        omp_set_lock(&next->lock);

        TreeNode* temp = curr;

        curr = next;

        next = (data <= curr->data) ? curr->left : curr->right;

        omp_unset_lock(&temp->lock);
    }

    TreeNode* new = createNode(data);
    curr->left  = (data <= curr->data) ? new : curr->left;
    curr->right = (data > curr->data) ? new : curr->right;

    omp_unset_lock(&curr->lock);

    return root;
}

// thought about also using searchNode but i would need to make some not pretty changes
TreeNode* deleteNode(TreeNode* root, int data) {
 
    if (root == NULL) return NULL;

    TreeNode* curr = root;
    TreeNode* next = root;

    // this is for coherent locking, sets next similar to other functions 
    omp_set_lock(&curr->lock);

    if (next != NULL && next->data != data) { next = (data <= curr->data) ? curr->left : curr->right; }
    // check and aquire lock if necessary
    while (next != NULL && next->data != data) {

        omp_set_lock(&next->lock);

        TreeNode* temp = curr;

        curr = next;

        next = (data <= curr->data) ? curr->left : curr->right;

        omp_unset_lock(&temp->lock);

        // node not found
        if (next == NULL) { omp_unset_lock(&curr->lock); return root; }
    }

    if (next != NULL && next->left == NULL && next->right == NULL) {

        curr->left  = (curr->left  != NULL && curr->left->data  == data) ? NULL : curr->left;
        curr->right = (curr->right != NULL && curr->right->data == data) ? NULL : curr->right;
        omp_unset_lock(&next->lock);
        // no node will come into undefined memory because we updated curr
        // omp_destroy_lock(&next->lock);
        free(next);
        if (next == root) return NULL;

    // if only one child
    } else if ((next != NULL && next->left != NULL && next->right == NULL) || (next != NULL && next->left == NULL && next->right != NULL)){

        TreeNode* child = (next->left != NULL) ? next->left : next->right;

        curr->left  = (curr->left != NULL  && curr->left->data == data)  ? child : curr->left;
        curr->right = (curr->right != NULL && curr->right->data == data) ? child : curr->right;
        omp_unset_lock(&next->lock);
        // no node will come into undefined memory because we updated curr
        // omp_destroy_lock(&next->lock);
        free(next);
        if (next == root) root = child;

    } else if (next != NULL) {

        TreeNode* succ = findMin(next->right);
        
        next->data = succ->data;
        // maybe we just deleted the child of next, we want updated tree
        next->right = deleteNode(next->right, succ->data);

    }

    omp_unset_lock(&curr->lock);

    return root;
}

bool searchNode(TreeNode* root, int data) {

    if (root == NULL) return false;
    if (root->data == data) return true;

    TreeNode* curr = root;
    omp_set_lock(&curr->lock);

    TreeNode* next = (data <= curr->data) ? curr->left : curr->right;

    while (next != NULL) {

        omp_set_lock(&next->lock);

        TreeNode* temp = curr;

        curr = next;

        if (curr->data == data) { omp_unset_lock(&curr->lock); omp_unset_lock(&temp->lock); return true; } 

        next = (data <= curr->data) ? curr->left : curr->right;

        omp_unset_lock(&temp->lock);
    }

    omp_unset_lock(&curr->lock);

    return false;
}

TreeNode* findMin(TreeNode* root) {

    if (root == NULL) return NULL;

    TreeNode* curr = root;
    omp_set_lock(&curr->lock);

    TreeNode* next = curr->left;

    while (next != NULL) {
        omp_set_lock(&next->lock);

        TreeNode* temp = curr;

        curr = next;

        next = curr->left;

        omp_unset_lock(&temp->lock);
    }

    omp_unset_lock(&curr->lock);

    return curr;
}

// there's no other way than locking all the tree if we want to print updated data
void inorderTraversal(TreeNode* root) {
    if (root == NULL) return;

    omp_set_lock(&root->lock);

    inorderTraversal(root->left);
    printf("%d\n", root->data);
    inorderTraversal(root->right);

    omp_unset_lock(&root->lock);
}

void preorderTraversal(TreeNode* root) {
    if (root == NULL) return;

    omp_set_lock(&root->lock);

    printf("%d\n", root->data);
    preorderTraversal(root->left);
    preorderTraversal(root->right);

    omp_unset_lock(&root->lock);
}

void postorderTraversal(TreeNode* root) {
    if (root == NULL) return;

    omp_set_lock(&root->lock);

    postorderTraversal(root->left);
    postorderTraversal(root->right);
    printf("%d\n", root->data);

    omp_unset_lock(&root->lock);
}

// assume it's not going to be called in parallel
void freeTree(TreeNode* root) {

    if (root == NULL) return;

    freeTree(root->left);
    freeTree(root->right);
    // omp_destroy_lock(&root->lock);
    free(root);
}
