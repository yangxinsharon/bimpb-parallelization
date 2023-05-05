/* Xin (Sharon) Yang
   SMU Mathematics
   Tree structure preconditioner */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "utilities.h"

#include "global_params.h"
#include "array.h"
#include "tree_node_struct.h"
#include "particle_struct.h"


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
extern double *h_pot;

extern double pi;
extern double one_over_4pi;
extern double bulk_coef;
extern double units_coef;
extern double epsw;
extern double epsp;
extern double eps;
extern double bulk_strength;  	
extern double kappa2;	
extern double kappa;
extern double **tr_xyz2D;

/* runtime treecode parameters */
static int s_numpars;
static int s_order;
static int s_max_per_leaf;
static double theta;

/* variables for tracking tree information */
static int s_min_level;
static int s_max_level;

/* global variables for reordering arrays */
static int *s_order_arr = NULL;
/* root node of tree */
static TreeNode *s_tree_root = NULL;



/* internal functions */
static int s_Setup(double *xyzminmax);
static int s_CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6],
                        int level);
static int s_PartitionEight(double xyzmms[6][8], double xl, double yl,
                            double zl, double lmax, double x_mid, double y_mid,
                            double z_mid, int ind[8][2]);



int TreecodeInitialization() {
	// transfer tr_xyz 1D to 2D
	int i,j;
	// double tr_xyz2D[3][nface];



	int s_numpars;
	s_numpars = nface;
	int s_min_level = 50000;
	int s_max_level = 0;
	int level = 0;

	double *xyzminmax;
	xyzminmax = (double*)calloc(6,sizeof(double));
	s_Setup(xyzminmax);
	s_tree_root = (TreeNode*)calloc(1, sizeof(TreeNode));
	s_CreateTree(s_tree_root, 0, s_numpars-1, xyzminmax, level);
	return 0;
}

/* preconditioning calculation */
int Nrow;
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

  	make_vector(xtemp, N);
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
  	free_vector(xtemp);
}


/**********************************************************/
int *psolve(double *z, double *r) {
/* r as original while z as scaled */
  
  	// clock_t start_p,finish_p;
  	// double total_p = 0;
  	// start_p = clock();

  	int i, j, idx = 0, nrow, nrow2, ibeg = 0, iend = 0;
  	int *ipiv, inc;
  	double **matrixA, *rhs;
  	double L1, L2, L3, L4, area;
  	double tp[3], tq[3], sp[3], sq[3];
  	double r_s[3], rs, irs, sumrs;
  	double G0, kappa_rs, exp_kappa_rs, Gk;
  	double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
  	double G10, G20, G1, G2, G3, G4;
  	double pre1, pre2;
	
  	pre1 = 0.5*(1.0+eps);
  	pre2 = 0.5*(1.0+1.0/eps);

  	/*
  	make_matrix(matrixA, 3, 3);
  	make_vector(rhs, 3);
  	matrixA[0][0] = 1.0; matrixA[0][1] = 2.0; matrixA[0][2] = 3.0;
  	matrixA[1][0] = 3.0; matrixA[1][1] = 2.0; matrixA[1][2] = 4.0;
  	matrixA[2][0] = 1.0; matrixA[2][1] = 5.0; matrixA[2][2] = 6.0;
  	rhs[0] = 1.0; rhs[1] = 2.0; rhs[2] = 3.0;
  	printf("%f %f %f\n",rhs[0],rhs[1],rhs[2]);
  	inc = lu_decomp( matrixA, 3, ipiv );
  	printf("%d %d %d\n",ipiv[0],ipiv[1],ipiv[2]);
  	lu_solve( matrixA, 3, ipiv, rhs );
  	printf("%f %f %f\n",rhs[0],rhs[1],rhs[2]);
  	free_matrix(matrixA);
  	free_vector(rhs);
  	exit(0);
  	*/

  	make_matrix(matrixA, 2*s_max_per_leaf, 2*s_max_per_leaf);
  	make_vector(ipiv, 2*s_max_per_leaf);
  	make_vector(rhs, 2*s_max_per_leaf);

  	while ( idx < s_numpars ) {
    	leaflength(s_tree_root, idx);
    	nrow  = Nrow;
    	nrow2 = nrow*2;
    	ibeg  = idx;
    	iend  = idx + nrow - 1;

    	for ( i = ibeg; i <= iend; i++ ) {
	      	tp[0] = tr_xyz[3*j]; //
	      	tp[1] = tr_xyz[3*j+1]; 
	      	tp[2] = tr_xyz[3*j+2]; 
      		tq[0] = tr_q[3*j]; //
      		tq[1] = tr_q[3*j+1]; //
      		tq[2] = tr_q[3*j+2]; //


      		for ( j = ibeg; j < i; j++ ) {
        		sp[0] = tr_xyz[3*j]; //s_particle_position[0][j];
        		sp[1] = tr_xyz[3*j+1]; //s_particle_position[1][j];
        		sp[2] = tr_xyz[3*j+2]; //s_particle_position[2][j];
        		sq[0] = tr_q[3*j]; //s_particle_normal[0][j];
        		sq[1] = tr_q[3*j+1]; //s_particle_normal[1][j];
        		sq[2] = tr_q[3*j+2]; //s_particle_normal[2][j];

        		r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
        		sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];

        		rs = sqrt(sumrs);
        		irs = 1.0/rs;
        		G0 = ONE_OVER_4PI * irs;
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
		
        		area = tr_area[j]; // s_particle_area[j];
		
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
	        	sp[0] = tr_xyz[3*j]; //s_particle_position[0][j];
	        	sp[1] = tr_xyz[3*j+1]; //s_particle_position[1][j];
	        	sp[2] = tr_xyz[3*j+2]; //s_particle_position[2][j];
	        	sq[0] = tr_q[3*j]; //s_particle_normal[0][j];
	        	sq[1] = tr_q[3*j+1]; //s_particle_normal[1][j];
	        	sq[2] = tr_q[3*j+2]; //s_particle_normal[2][j];
	        	
	        	r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
				sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
	        	rs = sqrt(sumrs);
	        	irs = 1.0/rs;
	        	G0 = ONE_OVER_4PI * irs;
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
      		rhs[i+nrow] = r[i+ibeg+s_numpars];
    	}

    	inc = lu_decomp( matrixA, nrow2, ipiv );
    	lu_solve( matrixA, nrow2, ipiv, rhs );

    	for ( i = 0; i < nrow; i++) {
      		z[i+ibeg] = rhs[i];
      		z[i+ibeg+s_numpars] = rhs[i+nrow];
    	}

    	//printf("%d %d %d %d\n", idx, ibeg, iend, nrow);

    	idx += nrow;

  	}
  	free_matrix(matrixA);
  	free_vector(rhs);
  	free_vector(ipiv);

  	// for ( i = 0; i < s_numpars; i++) {
  	//   z[i] = r[i]/pre1;
  	//   z[i+s_numpars] = r[i+s_numpars]/pre2;
  	// }
  	// finish_p = clock();
  	// total_p = (double)(finish_p - start_p);
  	// printf("psolve time is %f\n", total_p);
  	return 0;

}




int s_Setup(double *xyzminmax) {
/*	the smallest box containing the particles is determined. The particle
	postions and charges are copied so that they can be restored upon exit.
*/

	int i;
	// find bounds of Cartesian box enclosing the particles 
	xyzminmax[0] = MinVal(tr_xyz2D[0],s_numpars);
	xyzminmax[1] = MaxVal(tr_xyz2D[0],s_numpars);
	xyzminmax[2] = MinVal(tr_xyz2D[1],s_numpars);
	xyzminmax[3] = MaxVal(tr_xyz2D[1],s_numpars);
	xyzminmax[4] = MinVal(tr_xyz2D[2],s_numpars);
	xyzminmax[5] = MaxVal(tr_xyz2D[2],s_numpars);
	// xyzminmax[0] = MinVal(x,s_numpars);
	// xyzminmax[1] = MaxVal(x,s_numpars);
	// xyzminmax[2] = MinVal(y,s_numpars);
	// xyzminmax[3] = MaxVal(y,s_numpars);
	// xyzminmax[4] = MinVal(z,s_numpars);
	// xyzminmax[5] = MaxVal(z,s_numpars);

   // if ((orderarr=(int *) malloc(numpars*sizeof(int)))==NULL) {
	// 	printf("Error allocating copy variables!");
	// }
   	make_vector(s_order_arr, s_numpars);
	for (i=0; i<s_numpars; i++) {
		s_order_arr[i] = i;
	}
	return 0;
}



// Node P is input, which contains particles indexed from IBEG to IEND.
//fortran
// int create_tree(TreeNode *p, int ibeg, int iend, double *x, double *y, double *z, double *q, 
// 	int maxparnode, double xyzmm[6], int level, int numpars) {
// 	/*CREATE_TREE recursively create the tree structure. Node P is
// 	input, which contains particles indexed from IBEG to IEND. After
// 	the node parameters are set, subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
// 	Real array XYZMM contains the min and max values of the coordinates
// 	of the particle in P, thus defining the box.*/

// 	// local variables
// 	double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
// 	int ind[8][2]={0};
// 	double xyzmms[6][8];
// 	int i, j, limin, limax, err, loclev, numposchild;
// 	double lxyzmm[6];

// 	// allocate pointer
	
// 	// set node fields: number of particle, exist_ms and xyz bounds
// 	p->numpar = iend-ibeg+1;
// 	p->exist_ms = 0;
//    	p->x_min = xyzmm[0];
//    	p->x_max = xyzmm[1];
// 	p->y_min = xyzmm[2];
//    	p->y_max = xyzmm[3];
//    	p->z_min = xyzmm[4];
//    	p->z_max = xyzmm[5];
// 	// compute aspect ratio
//    	xl = p->x_max - p->x_min;
//    	yl = p->y_max - p->y_min;
//    	zl = p->z_max - p->z_min;


//     lmax = xl;
//     if (lmax < yl) lmax = yl;
//     if (lmax < zl) lmax = zl;

//     t1 = lmax;
//     t2 = xl;
//     if (t2 > yl) t2 = yl;
//     if (t2 > zl) t2 = zl;

//     if (t2 != 0.0) {
//         p->aspect = t1/t2;
//     } else {
//         p->aspect = 0.0;
//     }

// 	// midpoint coordinates , RADIUS and SQRADIUS
// 	p->x_mid = (p->x_max + p->x_min) / 2.0;
//    	p->y_mid = (p->y_max + p->y_min) / 2.0;
//    	p->z_mid = (p->z_max + p->z_min) / 2.0;
//    	t1 = p->x_max-p->x_mid;
//    	t2 = p->y_max-p->y_mid;
//    	t3 = p->z_max-p->z_mid;
//    	p->radius = sqrt(t1*t1+t2*t2+t3*t3);
// 	// set particle limits, tree level of node, and nullify children pointers
//    	p->ibeg = ibeg;
//    	p->iend = iend;
//    	p->level = level;

//    	if (s_max_level < level) {
//    	   s_max_level = level;
//    	}

//    	p->num_children = 0;

//    	make_vector(p->child, 8);
//    	for(i=0; i<8; i++) {
//    	    p->child[i] = (TreeNode*)calloc(1, sizeof(TreeNode));
//    	}


//    	if (p->numpar > s_max_per_leaf) {
// 		// set IND array to 0 and then call PARTITION routine.  IND array holds indices
// 		// of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
//    		xyzmms[0][0] = p->x_min;
//    		xyzmms[1][0] = p->x_max;
//    		xyzmms[2][0] = p->y_min;
//    		xyzmms[3][0] = p->y_max;
//    		xyzmms[4][0] = p->z_min;
//    		xyzmms[5][0] = p->z_max;

//    		for (i = 0; i < 8; i++) {
//             ind[i][0] = 0;
//             ind[i][1] = 0;
//         }

//    		ind[0][0] = ibeg;
//    		ind[0][1] = iend;
//    		x_mid = P->x_mid;
//    		y_mid = P->y_mid;
//    		z_mid = P->z_mid;
// 		partition_8(x,y,z,q,xyzmms,xl,yl,zl,lmax,&numposchild, x_mid,y_mid,z_mid,ind,numpars)	
// 	}

// 	// Shrink the box
// 	// if (1==2){
// 	// 	for (i=0; i<8; i++){
// 	// 		if (ind[i][0]<ind[i][1]) {
//     //             xyzmms[0][i] = minval(x[ind[i][0]:ind[i][1]])
//     //             xyzmms[1][i] = maxval(x[ind[i][0]:ind[i][1]])
//     //             xyzmms[2][i] = minval(y[ind[i][0]:ind[i][1]])
//     //             xyzmms[3][i] = maxval(y[ind[i][0]:ind[i][1]])
//     //             xyzmms[4][i] = minval(z[ind[i][0]:ind[i][1]])
//     //             xyzmms[5][i] = maxval(z[ind[i][0]:ind[i][1]])
// 	// 		}
// 	// 	}
// 	// }

// 	// create children if indicated and store info in parent
// 	loclev = level + 1;
// 	for (i=0; i<numposchild; i++) {
// 		if (ind[i][0] < ind[0][1]) {
// 			p->num_children = p->num_children + 1;
// 			for (j=0; i<6; j++){
// 				lxyzmm[j] = xyzmms[j][i];
// 			}
// 			create_tree(P->child[p->num_children-1]->p_to_tnode,ind[i][0], ind[i][1], 
// 				x, y, z, q, maxparnode, lxyzmm, loclev, numpars)
// 		}
// 	}


// }

/********************************************************/
int s_CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6], int level)
{
/*CREATE_TREE recursively create the tree structure. Node P is
  input, which contains particles indexed from IBEG to IEND. After
  the node parameters are set subdivision occurs if IEND-IBEG+1 > s_max_per_leaf.
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

    make_vector(p->child, 8);
    for (i = 0; i < 8; i++) {
        p->child[i] = (TreeNode*)calloc(1, sizeof(TreeNode));
    }


    if (p->numpar > s_max_per_leaf) {
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

        numposchild = s_PartitionEight(xyzmms, xl, yl, zl, lmax,
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
                s_CreateTree(p->child[p->num_children-1],
                             ind[i][0], ind[i][1], lxyzmm, loclev);
            }
        }
    } else {
        if (level < s_min_level) {
            s_min_level = level;
        }
    }

    return 0;
}
/**********************************************************/



// int partition_8(double *x,double *y,double *z,double *q,double *xyzmms[8],
// 	double xl,double yl,double zl,double lmax,int *numposchild,
// 	double x_mid, double y_mid,double z_mid,int *ind[2],int numpars) {
// 	/* PARTITION_8 determines the particle indices of the eight sub boxes 
// 	containing the particles after the box defined by particles I_BEG
// 	to I_END is divided by its midpoints in each coordinate direction. */
	
// 	int temp_ind, i;
// 	double critlen;

// 	numposchild = 1;
// 	critlen = lmax/sqrt(2.0);

// 	if (xl >= critlen) {
// 		temp_ind = partition(x,y,z,q,orderarr,ind[0][0],ind[0][1],x_mid);
//         ind[1][0] = temp_ind+1;
//         ind[1][1] = ind[0][1];
//         ind[0][1] = temp_ind;
//         for (i=0; i<6; i++) {
//         	xyzmms[i][1]=xyzmms[i][0];
//         }
//         xyzmms[1][0]=x_mid;
//         xyzmms[0][1]=x_mid;
//         numposchild=2*numposchild;
// 	}

// 	if (yl >= critlen) {
// 		for (i=0; i<numposchild; i++){

// 			partition(x,y,z,q,orderarr,ind[i][0],ind[i][1],y_mid,&temp_ind,numpars);
//         	ind[numposchild+i][0] = temp_ind+1;
//         	ind[numposchild+i][1] = ind[i][1];
//         	ind[i][1] = temp_ind;
//         	for (j=0; j<6;j++){
//         		xyzmms[j][numposchild+i] = xyzmms[j][i]
//         	}
//         	xyzmms[3][i] = y_mid;
//         	xyzmms[2][numposchild+i] = y_mid;
        	
//         	numposchild = 2*numposchild;
//         }
// 	}

// 	if (zl >= critlen) {
// 		for (i=0; i<numposchild; i++){
// 			partition(z,x,y,q,orderarr,ind[i][0],ind[i][1],z_mid,&temp_ind,numpars);
//             ind[numposchild+i][0] = temp_ind+1;
//             ind[numposchild+i][1] = ind[i][1];
//             ind[i][1]=temp_ind;
//             for (j=0; j<6;j++){
//             	xyzmms[j][numposchild+i] = xyzmms[j][i];
//             }
//             xyzmms[5][i] = z_mid;
//             xyzmms[4][numposchild+i] = z_mid;
// 		}
// 		numposchild = 2* numposchild;
// 	}


// }

/********************************************************/
static int s_PartitionEight(double xyzmms[6][8], double xl, double yl, 
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
        		lower = lower - 1;
        	}
        	if (upper != lower) {
        		a[upper] = a[lower];
            	b[upper] = b[lower];
            	c[upper] = c[lower];
            indarr[upper] = indarr[lower];
        	}
        	while ((upper < lower) && (val >= a[upper])) {	
        		upper = upper + 1;
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
   	}

   	else if (ibeg == iend) {
    	if (a[ibeg] <= val) {
    		midind = ibeg;
    	}
    	else {
    		midind = ibeg -1 ;
    	}
    }
   	else {
    	midind = ibeg -1;
   	}
   	return (midind); 
}
