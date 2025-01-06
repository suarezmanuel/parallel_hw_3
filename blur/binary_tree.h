// 340858935 Manuel Ignacio Suarez Del Solar
// 216270827 Azmon Avraham Zvi

#ifndef BINARY_TREE_H
#define BINARY_TREE_H

#include <stdbool.h>
#include <omp.h>

typedef struct TreeNode {
    int data;
    struct TreeNode *left;
    struct TreeNode *right;
    omp_lock_t lock;
} TreeNode;

TreeNode* createNode(int data);

TreeNode* insertNode(TreeNode* root, int data);

TreeNode* deleteNode(TreeNode* root, int data);

bool searchNode(TreeNode* root, int data);

TreeNode* findMin(TreeNode* root);

void inorderTraversal(TreeNode* root);

void preorderTraversal(TreeNode* root);

void postorderTraversal(TreeNode* root);

void freeTree(TreeNode* root);

#endif