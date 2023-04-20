/* Xin (Sharon) Yang
   SMU Mathematics
   Tree structure preconditioner */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gl_variables.h"
#include "gl_constants.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

// node structure and node type declarations
struct tnode; //incomplete declaration

typedef struct tnode_pointer {
	tnode *p_to_tnode;
}

typedef struct tnode {
	int numpar, ibeg, iend;
	double x_min, y_min, z_min;
	double x_max, y_max, z_max;
	double x_mid, y_mid, z_mid;
	double aspect, radius;
	int exist_ms, level, num_children;
	double ***ms; //dim(:,:,:,:)
	tnode_pointer child[8]; //dim(8)
};

tnode troot;


/* preconditioning calculation */
int Nrow;
void leaflength(tnode *P, int idx) {
	/* find the leaf length */

	if (idx == P->ibeg && P->num_children == 0 ) {
   	Nrow = p->numpar;
  	} 
  	else {
    	if ( p->num_children != 0 ) {
      	for ( int i = 0; i < p->num_children; i++ ){
        		leaflength(p->child[i], idx);
      	}
    	}
  	}
}



// /* set up Taylor expansion and its coefficients */
// void setup(double x[],double y[],double z[],double q[],int numpars,int order,int iflag,double *xyzminmax) {
// /*	SETUP allocates and initializes arrays needed for the Taylor expansion.
// 	Also, global variables are set and the Cartesian coordinates of
// 	the smallest box containing the particles is determined. The particle
// 	postions and charges are copied so that they can be restored upon exit.
// */

// 	// local variables
// 	int err, i, j, k;
// 	double t1;

// 	// global integers and reals:  TORDER, TORDERLIM and THETASQ
// 	int torder = order + 2;
// 	int torder2 = order;

// 	if (iflag != 1) {
// 		int orderoffset = 0;
// 	}
// 	else {
// 		int orderoffset = 1;
// 	}
// 	int torderlim = torder + orderoffset;

// 	// allocate global Taylor expansion variables
// 	if ((cf=(double *) malloc((torder+1)*sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	}

// 	if ((cf1=(double *) malloc(torderlim*sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	}

// 	if ((cf2=(double *) malloc(torderlim*sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	}


// 	if ((cf3=(double *) malloc(torderlim*sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	}

// 	if ((a=(double ***) calloc((torderlim+3),sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	} // a(-2:torderlim,-2:torderlim,-2:torderlim)

// 	if ((b=(double ***) calloc((torderlim+3),sizeof(double)))==NULL) {
// 		printf("Error allocationg Taylor variables!");
// 	} // b(-2:torderlim,-2:torderlim,-2:torderlim)


// 	for (i=0; i<torder+1; i++){
// 		cf[i] = double(i)+1.0;
// 	}
// 	for (i=1; i<torderlim+1; i++){
// 		t1 = 1.0/double(i);
// 		cf1[i] = t1;
// 		cf2[i] = 1.0-0.5*t1;
// 		cf3[i] = 1.0-t1;
// 	}

// 	// find bounds of Cartesian box enclosing the particles 
//     xyzminmax[0]=min(x[0:numpars]);
//     xyzminmax[1]=max(x[0:numpars]);
//     xyzminmax[2]=min(y[0:numpars]);
//     xyzminmax[3]=max(y[0:numpars]);
//     xyzminmax[4]=min(z[0:numpars]);
//     xyzminmax[5]=max(z[0:numpars]);

//     if ((orderarr=(int *) malloc(numpars*sizeof(int)))==NULL) {
// 		printf("Error allocating copy variables!");
// 	}

// 	for (i=1; i<numpars+1; i++) {
// 		orderarr[i] = i;
// 	}
// }



// Node P is input, which contains particles indexed from IBEG to IEND.
void create_tree(tnode *P, int ibeg, int iend, double *x, double *y, double *z, double *q, 
	int maxparnode, double xyzmm[6], int level, int numpars) {
	/*CREATE_TREE recursively create the tree structure. Node P is
	input, which contains particles indexed from IBEG to IEND. After
	the node parameters are set, subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
	Real array XYZMM contains the min and max values of the coordinates
	of the particle in P, thus defining the box.*/

	// local variables
	double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
	int ind[8][2]={0};
	double xyzmms[6][8];
	int i, j, limin, limax, err, loclev, numposchild;
	double lxyzmm[6];

	// allocate pointer
	
	// set node fields: number of particle, exist_ms and xyz bounds
	P->numpar = iend-ibeg+1;
	P->exist_ms = 0;
   P->x_min = xyzmm[0];
   P->x_max = xyzmm[1];
   P->y_min = xyzmm[2];
   P->y_max = xyzmm[3];
   P->z_min = xyzmm[4];
   P->z_max = xyzmm[5];
	/ compute aspect ratio
   xl = P->x_max - P->x_min;
   yl = P->y_max - P->y_min;
   zl = P->z_max - P->z_min;
   lmax = max(max(xl,yl),zl);
   t1 = lmax;
   t2 = min(min(xl,yl),zl);
   if (t2 != 0.0) { // t>1e-08?
      P->aspect = t1/t2;
   }
   else {
      P->aspect=0.0;
   }

	// midpoint coordinates , RADIUS and SQRADIUS
	P->x_mid = (P->x_max + P->x_min) / 2.0;
   P->y_mid = (P->y_max + P->y_min) / 2.0;
   P->z_mid = (P->z_max + P->z_min) / 2.0;
   t1 = P->x_max-P->x_mid;
   t2 = P->y_max-P->y_mid;
   t3 = P->z_max-P->z_mid;
   P->radius = sqrt(t1*t1+t2*t2+t3*t3);
	/ set particle limits, tree level of node, and nullify children pointers
   P->ibeg = ibeg;
   P->iend = iend;
   P->level = level;
   if (maxlevel < level) {
      maxlevel = level;
   }
   P->num_children = 0;
   for(i=0; i<8; i++) {
      P->child[i]->p_to_tnode = NULL;
   }
   if (P->numpar > maxparnode) {
	// set IND array to 0 and then call PARTITION routine.  IND array holds indices
	// of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
   		xyzmms[0][0] = P->x_min;
   		xyzmms[1][0] = P->x_max;
   		xyzmms[2][0] = P->y_min;
   		xyzmms[3][0] = P->y_max;
   		xyzmms[4][0] = P->z_min;
   		xyzmms[5][0] = P->z_max;

   		ind[0][0] = ibeg;
   		ind[0][1] = iend;
   		x_mid = P->x_mid;
   		y_mid = P->y_mid;
   		z_mid = P->z_mid;
		partition_8(x,y,z,q,xyzmms,xl,yl,zl,lmax,&numposchild, x_mid,y_mid,z_mid,ind,numpars)	
	}

	// Shrink the box
	// if (1==2){
	// 	for (i=0; i<8; i++){
	// 		if (ind[i][0]<ind[i][1]) {
    //             xyzmms[0][i] = minval(x[ind[i][0]:ind[i][1]])
    //             xyzmms[1][i] = maxval(x[ind[i][0]:ind[i][1]])
    //             xyzmms[2][i] = minval(y[ind[i][0]:ind[i][1]])
    //             xyzmms[3][i] = maxval(y[ind[i][0]:ind[i][1]])
    //             xyzmms[4][i] = minval(z[ind[i][0]:ind[i][1]])
    //             xyzmms[5][i] = maxval(z[ind[i][0]:ind[i][1]])
	// 		}
	// 	}
	// }

	// create children if indicated and store info in parent
	loclev = level + 1;
	for (i=0; i<numposchild; i++) {
		if (ind[i][0] < ind[0][1]) {
			p->num_children = p->num_children + 1;
			for (j=0; i<6; j++){
				lxyzmm[j] = xyzmms[j][i];
			}
			create_tree(P->child[p->num_children-1]->p_to_tnode,ind[i][0], ind[i][1], 
				x, y, z, q, maxparnode, lxyzmm, loclev, numpars)
		}
	}


}



void partition_8(double *x,double *y,double *z,double *q,double *xyzmms[8],
	double xl,double yl,double zl,double lmax,int *numposchild,
	double x_mid, double y_mid,double z_mid,int *ind[2],int numpars) {
	/* PARTITION_8 determines the particle indices of the eight sub boxes 
	containing the particles after the box defined by particles I_BEG
	to I_END is divided by its midpoints in each coordinate direction. */
	
	int temp_ind, i;
	double critlen;

	numposchild = 1;
	critlen = lmax/sqrt(2.0);

	if (xl >= critlen) {
		partition(x,y,z,q,orderarr,ind[0][0],ind[0][1],x_mid,&temp_ind,numpars);
        ind[1][0] = temp_ind+1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;
        for (i=0; i<6; i++) {
        	xyzmms[i][1]=xyzmms[i][0];
        }
        xyzmms[1][0]=x_mid;
        xyzmms[0][1]=x_mid;
        numposchild=2*numposchild;
	}

	if (yl >= critlen) {
		for (i=0; i<numposchild; i++){

			partition(x,y,z,q,orderarr,ind[i][0],ind[i][1],y_mid,&temp_ind,numpars);
        	ind[numposchild+i][0] = temp_ind+1;
        	ind[numposchild+i][1] = ind[i][1];
        	ind[i][1] = temp_ind;
        	for (j=0; j<6;j++){
        		xyzmms[j][numposchild+i] = xyzmms[j][i]
        	}
        	xyzmms[3][i] = y_mid;
        	xyzmms[2][numposchild+i] = y_mid;
        	
        	numposchild = 2*numposchild;
        }
	}

	if (zl >= critlen) {
		for (i=0; i<numposchild; i++){
			partition(z,x,y,q,orderarr,ind[i][0],ind[i][1],z_mid,&temp_ind,numpars);
            ind[numposchild+i][0] = temp_ind+1;
            ind[numposchild+i][1] = ind[i][1];
            ind[i][1]=temp_ind;
            for (j=0; j<6;j++){
            	xyzmms[j][numposchild+i] = xyzmms[j][i];
            }
            xyzmms[5][i] = z_mid;
            xyzmms[4][numposchild+i] = z_mid;
		}
		numposchild = 2* numposchild;
	}


}



void partition(double *a, double *b, double *c, double *q, int *indarr,
	int ibeg, int iend, double val, int *midind, int numpars) {
	/* PARTITION determines the index MIDIND, after partitioning
	in place the rrays A,B,C and Q, such that
	A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL.
	If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
	is returned as IBEG-1.*/

	// local variables 
	double ta, tb, tc, tq;
	int lower, upper, tind;

	if (ibeg < iend) {

		// temporarily store IBEG entries and set A(IBEG)=VAL for the partitoning algorithm.
		ta = a[ibeg];
        tb = b[ibeg];
        tc = c[ibeg];
        tq = q[ibeg];
        tind = indarr[ibeg];
        a[ibeg] = val;
        upper = ibeg;
        lower = iend;
        while (upper != lower) {
        	while ((upper < lower) && (val < a[lower])) {
        		lower = lower - 1;

        	}
        	if (upper != lower) {
        		a[upper] = a(lower);
               	b[upper] = b(lower);
               	c[upper] = c(lower);
               	q[upper] = q(lower);
               	indarr[upper] = indarr[lower];

        	}
        	while ((upper < lower) && (val >= a[upper])) {
        		upper = upper + 1;
        	}
        	if (upper != lower) {
        		a[lower] = a[upper];
               	b[lower] = b[upper];
               	c[lower] = c[upper];
               	q[lower] = q[upper];
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
        q[upper] = tq;
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
    
}




/* Prototypes */
int *matvec(double *alpha, double *x, double *beta, double *y);
void comp_soleng_wrapper(double soleng);
void comp_source_wrapper();
void comp_pot(const double* xvct, double *atmchr, double *chrpos, double *ptl, 
	double *tr_xyz,double *tr_q, double *tr_area, int nface, int nchr);
void comp_source( double* bvct, double *atmchr, double *chrpos, 
	double *tr_xyz, double *tr_q, int nface, int nchr);

 
void matvecmul(const double *x, double *y, double *q, int nface, 
	double *tr_xyz, double *tr_q, double *tr_area, double alpha, double beta) {
	int i, j;
	double pre1, pre2;
	double area, rs, irs, sumrs;
	double G0, kappa_rs, exp_kappa_rs, Gk;
	double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
	double G10, G20, G1, G2, G3, G4;
	double L1, L2, L3, L4;

    pre1=0.50*(1.0+eps); /* const eps=80.0 */
    pre2=0.50*(1.0+1.0/eps);
    for (i=0; i<nface; i++) {
    	double tp[3] = {tr_xyz[3*i], tr_xyz[3*i+1], tr_xyz[3*i+2]};
		double tq[3] = {tr_q[3*i], tr_q[3*i+1], tr_q[3*i+2]};

		double peng[2] = {0.0, 0.0};
		for (j=0; j<nface; j++) {
        	if (j != i) {
				double sp[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
				double sq[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
				double r_s[3] = {sp[0]-tp[0], sp[1]-tp[1], sp[2]-tp[2]};
				sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
				rs = sqrt(sumrs);
				irs = 1.0/sqrt(sumrs) ; //rsqrt(sumrs);
				G0 = one_over_4pi;
				G0 = G0*irs;
				kappa_rs = kappa*rs;
				exp_kappa_rs = exp(-kappa_rs);
				Gk = exp_kappa_rs*G0;
	
				cos_theta = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
				cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
	
				tp1 = G0*irs;
				tp2 = (1.0+kappa_rs) * exp_kappa_rs;
	
				G10 = cos_theta0*tp1;
				G20 = tp2*G10;
	
				G1 = cos_theta*tp1;
				G2 = tp2*G1;
	
				dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
				G3 = (dot_tqsq-3.0*cos_theta0*cos_theta) * irs*tp1;
				G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
				L1 = G1-eps*G2;							// K2
				L2 = G0-Gk;								// K1
				L3 = G4-G3;								// K4
				L4 = G10-G20/eps; //fdivide(G20,eps);	// K3
	
				double peng_old[2] = {x[j], x[j+nface]};
				area = tr_area[j];
				peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
				peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
        	}
		}

		y[i] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
		y[nface+i] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
	}

}



/* This subroutine wraps the matrix-vector multiplication */
int *matvec(double *alpha, double *x, double *beta, double *y) {
    matvecmul(x, y, tr_q, nface, tr_xyz, tr_q, tr_area, *alpha, *beta);
    return NULL;
}

/* This subroutine wraps the solvation energy computation */
void comp_soleng_wrapper(double soleng) {
    int i;
	double *chrptl;
	double units_para = 2.0;
    units_para = units_para *units_coef;
    units_para = units_para*pi;

	if ((chrptl=(double *) malloc(nface*sizeof(double)))==NULL) {
		printf("error in allcating chrptl");
	}

	comp_pot(xvct, atmchr, chrpos, chrptl, tr_xyz, tr_q, tr_area, nface, nchr);
	soleng=0.0;
	for (i=0; i<nface; i++) soleng = soleng+chrptl[i];
	soleng = soleng*units_para;
	printf("solvation energy = %f kcal/mol\n",soleng);
}



/* This subroutine calculates the element-wise potential */
void comp_pot(const double* xvct, double *atmchr, double *chrpos, double *ptl, 
	double *tr_xyz, double *tr_q, double *tr_area, int nface, int nchr) {
	int i, j;
    double sumrs, irs, rs, G0, Gk, kappa_rs, exp_kappa_rs;
    double cos_theta, G1, G2, L1, L2, tp1, tp2;
	for (j=0; j<nface; j++) {
    	ptl[j] = 0.0;
		double r[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
		double v[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
    	for (i=0; i<nchr; i++) {
        	double s[3] = {chrpos[3*i], chrpos[3*i+1], chrpos[3*i+2]};
			double r_s[3] = {r[0]-s[0], r[1]-s[1], r[2]-s[2]};
			sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
			rs = sqrt(sumrs);
			irs = 1.0/sqrt(sumrs);

        	G0 = one_over_4pi;
        	G0 = G0*irs;
        	kappa_rs = kappa*rs;
        	exp_kappa_rs = exp(-kappa_rs);
        	Gk = exp_kappa_rs*G0;

        	cos_theta = (v[0]*r_s[0]+v[1]*r_s[1]+v[2]*r_s[2]) * irs;

        	tp1 = G0*irs;
        	tp2 = (1.0+kappa_rs)*exp_kappa_rs;

        	G1 = cos_theta*tp1;
        	G2 = tp2*G1;

        	L1 = G1-eps*G2;
        	L2 = G0-Gk;

      		ptl[j] = ptl[j] + atmchr[i] * (L1*xvct[j]+L2*xvct[nface+j]) * tr_area[j];
		}
    }
}

/* This subroutine wraps the solvation energy computation */
void comp_source_wrapper() {
    comp_source(bvct, atmchr, chrpos, tr_xyz, tr_q, nface, nchr);
}


/* This subroutine calculates the source term of the integral equation */
/* atmchr=atom charge   chrpos=charge position */
/* bvct be located at readin.c */
void comp_source( double* bvct, double *atmchr, double *chrpos, 
	double *tr_xyz,double *tr_q, int nface, int nchr) {
	int i, j;
	double sumrs, cos_theta, irs, G0, G1, tp1;
	for (i=0; i<nface; i++) {
        bvct[i] = 0.0;
        bvct[i+nface] = 0.0;
        for (j=0; j<nchr; j++) {
            double r_s[3] = {chrpos[3*j]-tr_xyz[3*i], chrpos[3*j+1]-tr_xyz[3*i+1], 
            	chrpos[3*j+2]-tr_xyz[3*i+2]};
			sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2]; 
            cos_theta = tr_q[3*i]*r_s[0] + tr_q[3*i+1]*r_s[1] + tr_q[3*i+2]*r_s[2];
			irs = 1.0/sqrt(sumrs) ;//rsqrt(sumrs);//returns reciprocal square root of scalars and vectors.
            cos_theta = cos_theta*irs;
            G0 = one_over_4pi;//constant
            G0 = G0*irs;
            tp1 = G0*irs;
            G1 = cos_theta*tp1;
            bvct[i] = bvct[i]+atmchr[j]*G0;
            bvct[nface+i] = bvct[nface+i]+atmchr[j]*G1;
        }

    }
}
