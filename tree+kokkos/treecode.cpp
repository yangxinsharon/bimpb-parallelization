/* Xin (Sharon) Yang
   SMU Mathematics
   Tree structure preconditioner */

/* Inclusions */
/* c */
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
// #include <string.h>

/* c++ */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "utilities.h"
#include "array.h"
#include "tree_node_struct.h"
#include "gl_constants.h"

#include <Kokkos_Core.hpp>

extern int nface, nspt, natm, nchr;
extern int **extr_v;						//[3][nspt]
extern int **extr_f;						//[2][nface]
extern int **face, **face_copy;				//[3][nface]
extern double **vert, **snrm;				//[3][nspt];
extern double *tr_xyz, *tr_q;				//[3][nface]
extern double *tr_area, *bvct, *xvct;		//[nface];
extern double **atmpos;						//[3][natm/nchr];
extern double *atmrad, *atmchr, *chrpos;	//[natm/nchr]; 
extern double *work, *h;

extern double **tr_xyz2D, **tr_q2D;
extern int **leafarr;


// #ifdef __cplusplus
// extern "C"
// {
// #endif
// extern make_vector(v,n);

// extern double pi;
// extern double one_over_4pi;
// extern double bulk_coef;
// extern double units_coef;
// extern double epsw;
// extern double epsp;
// extern double eps;
// extern double bulk_strength;	//ion_strength in M
// extern double kappa2;	//kappa2=bulk_coef*bulk_strength/epsw;
// extern double kappa;

// extern int order;
// extern int maxparnode;
// extern double theta;

extern double** Make2DDoubleArray(int arraySizeX, int arraySizeY, char info[]);
extern int** Make2DIntArray(int arraySizeX, int arraySizeY,char info[]);
// #ifdef __cplusplus
// }
// #endif


/* variables for tracking tree information */
static int s_min_level;
static int s_max_level;

/* global variables used when computing potential/force */
static double s_target_position[3];
static double s_target_normal[3];

static double ***s_target_charge = NULL;
static double ***s_source_charge = NULL;

/* global variables for reordering arrays */
static int *s_order_arr = NULL;
/* root node of tree */
static TreeNode *s_tree_root = NULL;

static int Nleaf = 0;
static int Nleafc = 0;
static int Nrow;
/* internal functions */
int Setup(double xyz_limits[6]);
void leaflength(TreeNode *p, int idx);
int Partition(double *a, double *b, double *c, int *indarr,
	int ibeg, int iend, double val);
int CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6],int level);
int PartitionEight(double xyzmms[6][8], double xl, double yl,
	double zl, double lmax, double x_mid, double y_mid,
	double z_mid, int ind[8][2]);
int s_RunTreecode(TreeNode *p, double *tpoten_old,
                  double tempq[2][16], double peng[2]);
int s_ComputeTreePB(TreeNode *p, double tempq[2][16], double peng[2]);
int s_ComputeDirectPB(int ibeg, int iend, double *tpoten_old,
                             double peng[2]);
int RemoveNode(TreeNode *p);

/**********************************************************/
int TreecodeInitialization() {
    
    int level, i, j, k, mm, nn;

    /* variables needed for reorder */
    double *temp_area, *temp_source;
    double **temp_normal;

    double xyz_limits[6];
    printf("\nInitializing treecode...\n");

    /* setting variables global to file */
    // s_order = order;
    s_min_level = 50000;
    s_max_level = 0;

    level = 0;

    // make_matrix(temp_normal, 3, s_numpars);
    // make_vector(temp_area, nface);
    // make_vector(temp_source, 2 * nface);
    temp_normal = Make2DDoubleArray(3,nface,"temp_normal");
    temp_area=(double *) calloc(nface, sizeof(double));
    temp_source=(double *) calloc(2*nface, sizeof(double));


	Setup(xyz_limits);
	s_tree_root = (TreeNode*)calloc(1, sizeof(TreeNode));
	
	// s_CreateTree(s_tree_root, 0, nface-1, xyzminmax, level);
	CreateTree(s_tree_root, 0, nface-1, xyz_limits, level);
	printf ("Nleaf is %d\n",Nleaf);
    printf("Created tree for %d particles with max %d per node.\n\n",
           nface, maxparnode);
	
    memcpy(temp_normal[0], tr_q2D[0], nface*sizeof(double));
    memcpy(temp_normal[1], tr_q2D[1], nface*sizeof(double));
    memcpy(temp_normal[2], tr_q2D[2], nface*sizeof(double));
    memcpy(temp_area, tr_area, nface*sizeof(double));
    memcpy(temp_source, bvct, 2*nface*sizeof(double));

    for (i = 0; i < nface; i++) {
        tr_q2D[0][i]    	= temp_normal[0][s_order_arr[i]];
        tr_q2D[1][i]    	= temp_normal[1][s_order_arr[i]];
        tr_q2D[2][i]    	= temp_normal[2][s_order_arr[i]];
        tr_area[i]         	= temp_area[s_order_arr[i]];
        bvct[i]          	= temp_source[s_order_arr[i]];
        bvct[i + nface] = temp_source[s_order_arr[i] + nface];
    }

    // free_matrix(temp_normal);
    // free_vector(temp_area);
    // free_vector(temp_source);
	for(i=0;i<3;i++) {
		free(temp_normal[i]);
	}	
	free(temp_normal);
	free(temp_area);
	free(temp_source);

    // make_3array(s_target_charge, nface, 2, 16);
    // make_3array(s_source_charge, nface, 2, 16);

    // transform tr_xyz2D and tr_q2D to 1 dimension
	for (j=0; j<nface; j++){
		for (i=0; i<3; i++){
			tr_xyz[3*j+i] = tr_xyz2D[i][j];
			tr_q[3*j+i] = tr_q2D[i][j];
		}
	}

	return 0;
}

/********************************************************/
int TreecodeFinalization()
{

    int i;
    double *temp_area, *temp_source, *temp_xvct;
    double **temp_normal, **temp_position;

/***********reorder particles*************/

    // make_matrix(temp_position, 3, nface);
    // make_matrix(temp_normal, 3, nface);
    // make_vector(temp_area, nface);
    // make_vector(temp_source, 2 * nface);
    // make_vector(temp_xvct, 2 * nface);

    temp_position=Make2DDoubleArray(3,nface,"temp_position");
    temp_normal=Make2DDoubleArray(3,nface,"temp_normal");
    temp_area=(double *) calloc(nface, sizeof(double));
    temp_source=(double *) calloc(2*nface, sizeof(double));
    temp_xvct=(double *) calloc(2*nface, sizeof(double));

    memcpy(temp_position[0], tr_xyz2D[0], nface*sizeof(double));
    memcpy(temp_position[1], tr_xyz2D[1], nface*sizeof(double));
    memcpy(temp_position[2], tr_xyz2D[2], nface*sizeof(double));
    memcpy(temp_normal[0], tr_q2D[0], nface*sizeof(double));
    memcpy(temp_normal[1], tr_q2D[1], nface*sizeof(double));
    memcpy(temp_normal[2], tr_q2D[2], nface*sizeof(double));
    memcpy(temp_area, tr_area, nface*sizeof(double));
    memcpy(temp_source, bvct, 2*nface*sizeof(double));
    memcpy(temp_xvct, xvct, 2*nface*sizeof(double));

    for (i = 0; i < nface; i++) {
        tr_xyz2D[0][s_order_arr[i]]  = temp_position[0][i];
        tr_xyz2D[1][s_order_arr[i]]  = temp_position[1][i];
        tr_xyz2D[2][s_order_arr[i]]  = temp_position[2][i];
        tr_q2D[0][s_order_arr[i]]    = temp_normal[0][i];
        tr_q2D[1][s_order_arr[i]]    = temp_normal[1][i];
        tr_q2D[2][s_order_arr[i]]    = temp_normal[2][i];
        tr_area[s_order_arr[i]]      = temp_area[i];
        bvct[s_order_arr[i]] 		 = temp_source[i];
        bvct[s_order_arr[i] + nface] = temp_source[i + nface];
       	xvct[s_order_arr[i]]         = temp_xvct[i];
        xvct[s_order_arr[i] + nface] = temp_xvct[i + nface];
    }

    // free_matrix(temp_position);
    // free_matrix(temp_normal);
    // free_vector(temp_area);
    // free_vector(temp_source);
    // free_vector(temp_xvct);
    for(i=0;i<3;i++) {
		free(temp_position[i]);
	}	
	free(temp_position);
    for(i=0;i<3;i++) {
    	free(temp_normal[i]);
    }
	free(temp_normal);

    free(temp_area);
    free(temp_source);
    free(temp_xvct);


/***********clean tree structure**********/
    // free_3array(s_target_charge);
    // free_3array(s_source_charge);

    RemoveNode(s_tree_root);
    free(s_tree_root);

    free_vector(s_order_arr);
/*****************************************/

    printf("\nTABIPB tree structure has been deallocated.\n\n");

    return 0;
}

/********************************************************/
int RemoveNode(TreeNode *p)
{
/* REMOVE_NODE recursively removes each node from the
 * tree and deallocates its memory for MS array if it exits. */
    int i;

    if (p->num_children > 0) {
        for (i = 0; i < 8; i++) {
            RemoveNode(p->child[i]);
            free(p->child[i]);
        }
        free(p->child);
    }

    return 0;
}


/********************************************************/
/* preconditioning calculation */
// int Nrow;
void leaflength(TreeNode *p, int idx) {
	/* find the leaf length */
	int i;
	if (idx == p->ibeg && p->num_children == 0 ) {
   		Nrow = p->numpar;
  	} else {
    	if ( p->num_children != 0 ) {
      		for ( i = 0; i < p->num_children; i++ ){
        		leaflength(p->child[i], idx);
      		}
    	}
  	}

}


/**********************************************************/
/* lapack provide lu decomposition, however, something    */
/* is wrong with cmake ************************************/
/**********************************************************/
int lu_decomp( double **A, int N, int *ipiv ) {

	int i, j, k, imax;
	double maxA, *ptr, absA, Tol = 1.0e-14;

  	for ( i = 0; i <= N; i++ ){
   		ipiv[i] = i; // record pivoting number
  	}

  	for ( i = 0; i < N; i++ ) {
   		maxA = 0.0;
   		imax = i;
   		for (k = i; k < N; k++){
   	  		if ((absA = fabs(A[k][i])) > maxA) {
   	   			maxA = absA;
   	    		imax = k;
   			}
   	  	}

   		if (maxA < Tol) return 0; //failure, matrix is degenerate	
   		
   		if (imax != i) {
	   	  	//pivoting P
	   	  	j = ipiv[i];
	   	  	ipiv[i] = ipiv[imax];
	   	  	ipiv[imax] = j;	
	   	  	//pivoting rows of A
	   	  	ptr = A[i];
	   	  	A[i] = A[imax];
	   	  	A[imax] = ptr;	
	   	  	//counting pivots starting from N (for determinant)
	   	  	ipiv[N]++;
	   	}	
   		for (j = i + 1; j < N; j++) {
   	  		A[j][i] /= A[i][i];	
   	  		for (k = i + 1; k < N; k++){
   	  	 		A[j][k] -= A[j][i] * A[i][k];
   	  		}
   		}
  	}

  	return 1;
}

void lu_solve( double **matrixA, int N, int *ipiv, double *rhs ) {
  	/* b will contain the solution */
  	double *xtemp;

  	// make_vector(xtemp, N);
  	xtemp=(double *) calloc(N, sizeof(double));
  	int i, k ;
  	for (i = 0; i < N; i++) {
   		xtemp[i] = rhs[ipiv[i]];

   		for (k = 0; k < i; k++){
      		xtemp[i] -= matrixA[i][k] * xtemp[k];
   		}
  	}

  	for (i = N - 1; i >= 0; i--) {
    	for (k = i + 1; k < N; k++){
      		xtemp[i] -= matrixA[i][k] * xtemp[k];
    	}

    	xtemp[i] = xtemp[i] / matrixA[i][i];
  	}

  	for (i = 0; i < N; i++) {
    	rhs[i] = xtemp[i];
  	}
  	// free_vector(xtemp);
  	free(xtemp);
}


/**********************************************************/
int *psolve(double *z, double *r) {
/* r as original while z as scaled */

  	int i, j, idx = 0, nrow, nrow2, ibeg = 0, iend = 0;
  	int *ipiv, inc;
  	double **matrixA; 
  	double *rhs;
  	double L1, L2, L3, L4, area;
  	double tp[3], tq[3], sp[3], sq[3];
  	double r_s[3], rs, irs, sumrs;
  	double G0, kappa_rs, exp_kappa_rs, Gk;
  	double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
  	double G10, G20, G1, G2, G3, G4;
  	double pre1, pre2;
	
  	pre1 = 0.5*(1.0+eps);
  	pre2 = 0.5*(1.0+1.0/eps);

  	// make_matrix(matrixA, 2*maxparnode, 2*maxparnode);
  	// make_vector(ipiv, 2*maxparnode);
  	// make_vector(rhs, 2*maxparnode);
  	matrixA=Make2DDoubleArray(2*maxparnode, 2*maxparnode, "matrixA");
	ipiv=(int *) calloc(2*maxparnode, sizeof(int));
	rhs=(double *) calloc(2*maxparnode, sizeof(double));
  	printf("maxparnode is %d\n", maxparnode);

//

	int arridx = 0, inleaf = 0;
	// inleaf = ceil(nface/maxparnode);
	// printf("inleaf is %d",inleaf);
	int leafarr[3][Nleaf];

	// leafarr=Make2DIntArray(3,Nleaf,"leafarr");


	while ( idx < nface ) {
	    leaflength(s_tree_root, idx);

	    nrow  = Nrow;
	    nrow2 = nrow*2;
	    ibeg  = idx;
	    iend  = idx + nrow - 1;
	    leafarr[0][arridx] = ibeg;
	    // printf("ibeg is: %d\n",leafarr[0][arridx]);

	    leafarr[1][arridx] = iend;
	    printf("iend is: %d\n",leafarr[1][arridx]);
	    leafarr[2][arridx] = nrow;
	    // printf("nrow is: %d\n",leafarr[2][arridx]);
	    arridx += 1;
	    Nleafc += 1;
	    // printf("arridx is: %d\n",arridx);
	    // printf("idx ibeg iend nrow: %d, %d, %d, %d\n",idx, leafarr[0][arridx], leafarr[1][arridx],leafarr[2][arridx] );
		// printf("idx ibeg iend nrow is %d, %d, %d, %d\n",idx,ibeg,iend,nrow);
		idx += nrow;
	}

	system("pause");


  	while ( idx < nface ) {
    	leaflength(s_tree_root, idx);

    	nrow  = Nrow;
    	nrow2 = nrow*2;
    	ibeg  = idx;
    	iend  = idx + nrow - 1;
    	Nleafc += 1;
    	// printf("idx ibeg iend is %d, %d, %d\n",idx,ibeg,iend);


    	for ( i = ibeg; i <= iend; i++ ) {
    		tp[0] = tr_xyz2D[0][i];
			tp[1] = tr_xyz2D[1][i];
			tp[2] = tr_xyz2D[2][i];
			tq[0] = tr_q2D[0][i];
			tq[1] = tr_q2D[1][i];
			tq[2] = tr_q2D[2][i];

      		for ( j = ibeg; j < i; j++ ) {
        		sp[0] = tr_xyz2D[0][j];
        		sp[1] = tr_xyz2D[1][j];
        		sp[2] = tr_xyz2D[2][j];
        		sq[0] = tr_q2D[0][j];
        		sq[1] = tr_q2D[1][j];
        		sq[2] = tr_q2D[2][j];    			
				
        		r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
        		sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];

        		rs = sqrt(sumrs);
        		irs = 1.0/rs;
        		G0 = one_over_4pi * irs;
        		kappa_rs = kappa * rs; //
        		exp_kappa_rs = exp(-kappa_rs);
        		Gk = exp_kappa_rs * G0;
		
        		cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
        		cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
        		tp1 = G0* irs;
        		tp2 = (1.0 + kappa_rs) * exp_kappa_rs;
		
        		G10 = cos_theta0 * tp1;
        		G20 = tp2 * G10;
		
        		G1 = cos_theta * tp1;
        		G2 = tp2 * G1;
		
        		dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
        		G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
        		G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
		
        		area = tr_area[j]; 
		
        		L1 = G1 - eps*G2;
        		L2 = G0 - Gk;
        		L3 = G4 - G3;
        		L4 = G10 - G20/eps;

        		matrixA[i-ibeg][j-ibeg] = -L1*area;
        		matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
        		matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
        		matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
      		}

      		matrixA[i-ibeg][i-ibeg] = pre1;
      		matrixA[i+nrow-ibeg][i+nrow-ibeg] = pre2;

      		for ( j = i+1; j <= iend; j++ ) {
        		sp[0] = tr_xyz2D[0][j];
        		sp[1] = tr_xyz2D[1][j];
        		sp[2] = tr_xyz2D[2][j];
        		sq[0] = tr_q2D[0][j];
        		sq[1] = tr_q2D[1][j];
        		sq[2] = tr_q2D[2][j];      			

	        	r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
				sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
	        	rs = sqrt(sumrs);
	        	irs = 1.0/rs;
	        	G0 = one_over_4pi * irs;
	        	kappa_rs = kappa * rs;
	        	exp_kappa_rs = exp(-kappa_rs);
	        	Gk = exp_kappa_rs * G0;
	        	
	        	cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
	        	cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
	        	tp1 = G0* irs;
	        	tp2 = (1.0 + kappa_rs) * exp_kappa_rs;
		
	        	G10 = cos_theta0 * tp1;
	        	G20 = tp2 * G10;
			
	        	G1 = cos_theta * tp1;
	        	G2 = tp2 * G1;
			
	        	dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
	        	G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
	        	G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
	        	area = tr_area[j];
		
	        	L1 = G1 - eps*G2;
	        	L2 = G0 - Gk;
	        	L3 = G4 - G3;
	        	L4 = G10 - G20/eps;
		
	        	matrixA[i-ibeg][j-ibeg] = -L1*area;
	        	matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
	        	matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
	        	matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
      		}
    	}

    	for ( i = 0; i < nrow; i++) {
      		rhs[i] = r[i+ibeg];
      		rhs[i+nrow] = r[i+ibeg+nface];
    	}

    	inc = lu_decomp( matrixA, nrow2, ipiv );
    	lu_solve( matrixA, nrow2, ipiv, rhs );

    	for ( i = 0; i < nrow; i++) {
      		z[i+ibeg] = rhs[i];
      		z[i+ibeg+nface] = rhs[i+nrow];
    	}

    	//printf("%d %d %d %d\n", idx, ibeg, iend, nrow);

    	idx += nrow;

  	}
  	printf("Nleafc is %d\n",Nleafc);
  	// free_matrix(matrixA);
  	// free_vector(rhs);
  	// free_vector(ipiv);

    for(i=0;i<2*maxparnode;i++) {
		free(matrixA[i]);
	}	
	free(matrixA);

  	free(rhs);
  	free(ipiv);

  	// for ( i = 0; i < nface; i++) {
  	//   z[i] = r[i]/pre1;
  	//   z[i+nface] = r[i+nface]/pre2;
  	// }

  	return 0;

}



/**********************************************************/
// int s_Setup(double *xyzminmax) {
int Setup(double xyz_limits[6]) {
/*	the smallest box containing the particles is determined. The particle
	postions and charges are copied so that they can be restored upon exit.
*/

	int i;
	// find bounds of Cartesian box enclosing the particles 
	xyz_limits[0] = MinVal(tr_xyz2D[0],nface);
	xyz_limits[1] = MaxVal(tr_xyz2D[0],nface);
	xyz_limits[2] = MinVal(tr_xyz2D[1],nface);
	xyz_limits[3] = MaxVal(tr_xyz2D[1],nface);
	xyz_limits[4] = MinVal(tr_xyz2D[2],nface);
	xyz_limits[5] = MaxVal(tr_xyz2D[2],nface);

   // if ((orderarr=(int *) malloc(numpars*sizeof(int)))==NULL) {
	// 	printf("Error allocating copy variables!");
	// }
   	// make_vector(s_order_arr, nface);
   	s_order_arr=(int *) calloc(nface, sizeof(int));

	for (i=0; i<nface; i++) {
		s_order_arr[i] = i;
	}
	return 0;
}


/********************************************************/
int CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6], int level)
{
/*CREATE_TREE recursively create the tree structure. Node P is
  input, which contains particles indexed from IBEG to IEND. After
  the node parameters are set subdivision occurs if IEND-IBEG+1 > maxparnode.
  Real array XYZMM contains the min and max values of the coordinates
  of the particle in P, thus defining the box. */

  /* local variables */
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int ind[8][2];
    double xyzmms[6][8];
    int i, j, loclev, numposchild;
    double lxyzmm[6];

/* set node fields: number of particles, exist_ms and xyz bounds */
    p->numpar = iend-ibeg+1;
    p->exist_ms = 0;

    p->x_min = xyzmm[0];
    p->x_max = xyzmm[1];
    p->y_min = xyzmm[2];
    p->y_max = xyzmm[3];
    p->z_min = xyzmm[4];
    p->z_max = xyzmm[5];

/* compute aspect ratio */
    xl = p->x_max-p->x_min;
    yl = p->y_max-p->y_min;
    zl = p->z_max-p->z_min;

    lmax = xl;
    if (lmax < yl) lmax = yl;
    if (lmax < zl) lmax = zl;

    t1 = lmax;
    t2 = xl;
    if (t2 > yl) t2 = yl;
    if (t2 > zl) t2 = zl;

    if (t2 != 0.0) {
        p->aspect = t1/t2;
    } else {
        p->aspect = 0.0;
    }

/* midpoint coordinates, RADIUS and SQRADIUS */
    p->x_mid = (p->x_max + p->x_min) / 2.0;
    p->y_mid = (p->y_max + p->y_min) / 2.0;
    p->z_mid = (p->z_max + p->z_min) / 2.0;
    t1 = p->x_max - p->x_mid;
    t2 = p->y_max - p->y_mid;
    t3 = p->z_max - p->z_mid;
    p->radius = sqrt(t1*t1 + t2*t2 + t3*t3);

/* set particle limits, tree level of node, and nullify children pointers */
    p->ibeg = ibeg;
    p->iend = iend;
    p->level = level;

    if (s_max_level < level) s_max_level = level;

    p->num_children = 0;

    // make_vector(p->child, 8);
    p->child=(sTreeNode **) calloc(8, sizeof(sTreeNode*));
    for (i = 0; i < 8; i++) {
        p->child[i] = (TreeNode*)calloc(1, sizeof(TreeNode));
    }


    if (p->numpar > maxparnode) {
/* set IND array to 0 and then call PARTITION routine. IND array holds indices
 * of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1 */

        xyzmms[0][0] = p->x_min;
        xyzmms[1][0] = p->x_max;
        xyzmms[2][0] = p->y_min;
        xyzmms[3][0] = p->y_max;
        xyzmms[4][0] = p->z_min;
        xyzmms[5][0] = p->z_max;

        for (i = 0; i < 8; i++) {
            ind[i][0] = 0;
            ind[i][1] = 0;
        }

        ind[0][0] = ibeg;
        ind[0][1] = iend;
        x_mid = p->x_mid;
        y_mid = p->y_mid;
        z_mid = p->z_mid;

        numposchild = PartitionEight(xyzmms, xl, yl, zl, lmax,
        	x_mid, y_mid, z_mid, ind);

/* Shrink the box */
        for (i = 0; i < 8; i++) {
            if (ind[i][0] < ind[i][1]) {
                xyzmms[0][i] = MinVal(&tr_xyz2D[0][ind[i][0]], ind[i][1]-ind[i][0]);
                xyzmms[1][i] = MaxVal(&tr_xyz2D[0][ind[i][0]], ind[i][1]-ind[i][0]);
                xyzmms[2][i] = MinVal(&tr_xyz2D[1][ind[i][0]], ind[i][1]-ind[i][0]);
                xyzmms[3][i] = MaxVal(&tr_xyz2D[1][ind[i][0]], ind[i][1]-ind[i][0]);
                xyzmms[4][i] = MinVal(&tr_xyz2D[2][ind[i][0]], ind[i][1]-ind[i][0]);
                xyzmms[5][i] = MaxVal(&tr_xyz2D[2][ind[i][0]], ind[i][1]-ind[i][0]);
            }
        }
/* create children if indicated and store info in parent */
        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                p->num_children = p->num_children + 1;
                for (j = 0; j < 6; j++) {
                    lxyzmm[j] = xyzmms[j][i];
                }
                CreateTree(p->child[p->num_children-1],
                	ind[i][0], ind[i][1], lxyzmm, loclev);
            }
        }
    } 
    else {
    	Nleaf += 1;
        if (level < s_min_level) {
            s_min_level = level;
        }
    }

    // leaf_arr = Make2DIntArray(2,Nleaf,)
    
    return 0;
}

/********************************************************/
// int RunTreecode(TreeNode *p, double *tpoten_old, double tempq[2][16],double peng[2])
// {
//   /* RunTreecode() is self recurrence function */
//     double tx, ty, tz, dist, pengchild[2];
//     int i;


//   /* determine DISTSQ for MAC test */
//     tx = p->x_mid - s_target_position[0];
//     ty = p->y_mid - s_target_position[1];
//     tz = p->z_mid - s_target_position[2];
//     dist = sqrt(tx*tx + ty*ty + tz*tz);

//   /* initialize potential energy */
//     peng[0] = 0.0;
//     peng[1] = 0.0;

// /* If MAC is accepted and there is more than 1 particale in the */
// /* box use the expansion for the approximation. */

//     if (p->radius < dist*theta && p->numpar > 40) {
//         s_ComputeTreePB(p, tempq, peng);
//     } else {
//         if (p->num_children == 0) {
//             s_ComputeDirectPB(p->ibeg, p->iend, tpoten_old, peng);
//         } else {
//       /* If MAC fails check to see if there are children. If not, perform */
//       /* direct calculation.  If there are children, call routine */
//       /* recursively for each. */
//             for (i = 0; i < p->num_children; i++) {
//                 pengchild[0] = 0.0;
//                 pengchild[1] = 0.0;
//                 s_RunTreecode(p->child[i], tpoten_old, tempq, pengchild);
//                 peng[0] += pengchild[0];
//                 peng[1] += pengchild[1];
//             }
//         }
//     }

//     return 0;
// }

/********************************************************/
int PartitionEight(double xyzmms[6][8], double xl, double yl, 
	double zl, double lmax, double x_mid, double y_mid, double z_mid, int ind[8][2]) {
/* PARTITION_8 determines the particle indices of the eight sub boxes
 * containing the particles after the box defined by particles I_BEG
 * to I_END is divided by its midpoints in each coordinate direction.
 * The determination of the indices is accomplished by the subroutine
 * PARTITION. A box is divided in a coordinate direction as long as the
 * resulting aspect ratio is not too large. This avoids the creation of
 * "narrow" boxes in which Talyor expansions may become inefficient.
 * On exit the INTEGER array IND (dimension 8 x 2) contains
 * the indice limits of each new box (node) and NUMPOSCHILD the number
 * of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
 * that box J is empty.*/
    int temp_ind, i, j;
    double critlen;
    int numposchild;

    numposchild = 1;
    critlen = lmax/sqrt(2.0);


    if (xl >= critlen) {
        temp_ind = Partition(tr_xyz2D[0],tr_xyz2D[1],tr_xyz2D[2], s_order_arr, ind[0][0], ind[0][1], x_mid);
        ind[1][0] = temp_ind+1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;
        for (i = 0; i < 6; i++) {
            xyzmms[i][1] = xyzmms[i][0];
        }
        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        numposchild *= 2;
    }

    if (yl >= critlen) {
        for (i = 0; i < numposchild; i++) {
            temp_ind = Partition(tr_xyz2D[1],tr_xyz2D[0],tr_xyz2D[2], s_order_arr, ind[i][0], ind[i][1], y_mid);
            ind[numposchild+i][0] = temp_ind+1;
            ind[numposchild+i][1] = ind[i][1];
            ind[i][1] = temp_ind;
            for (j = 0; j < 6; j++) {
                xyzmms[j][numposchild+i] = xyzmms[j][i];
            }
            xyzmms[3][i] = y_mid;
            xyzmms[2][numposchild+i] = y_mid;
        }
        numposchild *= 2;
    }

    if (zl >= critlen) {
        for (i = 0; i < numposchild; i++) {
            temp_ind = Partition(tr_xyz2D[2],tr_xyz2D[0],tr_xyz2D[1], s_order_arr, ind[i][0], ind[i][1], z_mid);
            ind[numposchild+i][0] = temp_ind+1;
            ind[numposchild+i][1] = ind[i][1];
            ind[i][1] = temp_ind;
            for (j = 0; j < 6; j++) {
                xyzmms[j][numposchild+i] = xyzmms[j][i];
            }
            xyzmms[5][i] = z_mid;
            xyzmms[4][numposchild+i] = z_mid;
        }
        numposchild *= 2;
    }

    return (numposchild);
}

/********************************************************/
int Partition(double *a, double *b, double *c, int *indarr,
	int ibeg, int iend, double val) {
	/* PARTITION determines the index MIDIND, after partitioning
	in place the arrays A,B,C and Q, such that
	A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL.
	If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
	is returned as IBEG-1.*/

	// local variables 
	double ta, tb, tc;
	int lower, upper, tind;
	int midind;

	if (ibeg < iend) {

		// temporarily store IBEG entries and set A(IBEG)=VAL for the partitoning algorithm.
		ta = a[ibeg];
      	tb = b[ibeg];
      	tc = c[ibeg];
      	tind = indarr[ibeg];
      	a[ibeg] = val;
      	upper = ibeg;
      	lower = iend;

      	while (upper != lower) {
      		while ((upper < lower) && (val < a[lower])) {
        		lower -= 1;
        	}
        	if (upper != lower) {
        		a[upper] = a[lower];
            	b[upper] = b[lower];
            	c[upper] = c[lower];
            	indarr[upper] = indarr[lower];
        	}
        	while ((upper < lower) && (val >= a[upper])) {	
        		upper += 1;
        	}
        	if (upper != lower) {
        		a[lower] = a[upper];
         		b[lower] = b[upper];
            	c[lower] = c[upper];
            	indarr[lower] = indarr[upper];
        	}
    	}
     	midind = upper;

		// replace TA in position UPPER and change MIDIND if TA > VAL
      	if (ta > val) {
      		midind = upper - 1;
      	}
      	a[upper] = ta;
      	b[upper] = tb;
      	c[upper] = tc;
      	indarr[upper] = tind;
   	} else if (ibeg == iend) {
    	if (a[ibeg] <= val) {
    		midind = ibeg;
    	} else {
    		midind = ibeg -1 ;
    	}
    } else {
    	midind = ibeg -1;
   	}
   	return (midind); 
}
