/**************************************************************************
* FILE NAME: tree_node_struct.h                                           *
*                                                                         *
* PURPOSE: typedef of stuct used to store node of treecode                *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 01/12/2018  Leighton Wilson   Created, moved from treecode header       *
*                                                                         *
**************************************************************************/

#ifndef H_TREE_NODE_STRUCT_H
#define H_TREE_NODE_STRUCT_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct sTreeNode {

    int node_idx;
    int numpar, ibeg, iend;
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    double radius, aspect;
    int level, num_children, exist_ms;
    double **ms;
    struct sTreeNode **child;

} TreeNode;


#ifdef __cplusplus
}
#endif

#endif /* H_TREE_NODE_STRUCT */

// yang:
// // node structure and node type declarations
// struct tnode; //incomplete declaration

// typedef struct tnode_pointer {
//     tnode *p_to_tnode;
// }

// typedef struct tnode {
//     int numpar, ibeg, iend;
//     double x_min, y_min, z_min;
//     double x_max, y_max, z_max;
//     double x_mid, y_mid, z_mid;
//     double aspect, radius;
//     int exist_ms, level, num_children;
//     double ***ms; //dim(:,:,:,:)
//     tnode_pointer child[8]; //dim(8)
// };

// tnode *troot = NULL;


