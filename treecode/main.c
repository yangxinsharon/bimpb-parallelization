/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pp_timer.h"


#include "array.h"
#include <time.h> // yang
#include "gl_constants.h"

extern int nface, nspt, natm, nchr;			// number of faces, points, atoms, and charges
extern int **extr_v;								// [3][nspt]
extern int **extr_f;								// [2][nface]
extern int **face, **face_copy;				// [3][nface]
extern double **vert, **snrm;					// [3][nspt]
extern double *tr_xyz, *tr_q;					// [3][nface]
extern double *tr_area, *bvct, *xvct;		// [nface]
extern double **atmpos;							// [3][natm/nchr]
extern double *atmrad, *atmchr, *chrpos;	// [natm/nchr]
extern double *work, *h;
extern double *h_pot;
extern const double eps;
extern tnode troot;

int main(int argc, char *argv[]) {
	/*variables local to main*/
	int i,j;
	double s[3], pot=0.0, sum=0.0, pot_temp=0.0;	// potential
	double ptl, soleng, t1, t2;
	char fname[16], density[16];
	extern void readin(char fname[16], char density[16]);

	/*GMRES related variables*/
	static long int info;
	long int RESTRT, ldw, ldh, iter, N;
	double resid;
	extern void comp_source_wrapper();				
	extern void comp_soleng_wrapper(double soleng);	
	extern int *matvec(double *alpha, double *x, double *beta, double *y); 
	extern int *psolve(double *z, double *r); 
	extern int gmres_(long int *n, double *b, double *x, long int *restrt, double *work, long int *ldw, 
		double *h, long int *ldh, long int *iter, double *resid, int *matvec (), int *psolve (), long int *info);

   extern void timer_start(char *n); 
   extern void timer_end(void); 
	timer_start("TOTAL_TIME");
	printf("%d %s %s %s \n", argc, argv[0], argv[1], argv[2]);

	/* read in structural information */
   // sprintf(fname, "1ajj");
   // sprintf(density, "1");
   sprintf(fname,"%s",argv[1]);
   sprintf(density,"%s",argv[2]);
	readin(fname, density);
	comp_source_wrapper(); //wraps the solvation energy computation

	/* parameters for GMRES */
	RESTRT=10;
	N=2*nface;
	ldw=N;
	ldh=RESTRT+1;
	iter=100;
	resid=1e-4;
	xvct=(double *) calloc(N, sizeof(double));

	work=(double *) calloc (ldw*(RESTRT+4), sizeof(double));
	h=(double *) calloc (ldh*(RESTRT+2), sizeof(double));

	gmres_(&N, bvct, xvct, &RESTRT, work, &ldw, h, &ldh, &iter, &resid, &matvec, &psolve, &info);

	soleng=0.0;

	comp_soleng_wrapper(soleng); //wraps the solvation energy computation
	timer_end();

	/* free memory */
	for(i=0;i<3;i++) {
		free(extr_v[i]);
	}
	free(extr_v);

	for(i=0;i<3;i++) {
		free(vert[i]);
	}
	free(vert);

	for(i=0;i<3;i++) {
		free(snrm[i]);
	}
	free(snrm);

	for(i=0;i<3;i++) {
		free(face[i]);
	}
	free(face);

	for(i=0;i<3;i++) {
		free(extr_f[i]);
	}
	free(extr_f);

	for(i=0;i<3;i++) {
		free(atmpos[i]);
	}
	free(atmpos);

	free(tr_xyz);
	free(tr_q);

	free(tr_area);
	free(bvct);
	free(xvct);
	free(atmchr);
	free(chrpos);

   return 0;
}

// ****************************************************************
// int *psolve(double *z, double *r) {
// /*r as original while z as scaled*/
// 	int i;
// 	double scale1,scale2;
// 	scale1=0.5*(1.0+eps);
// 	scale2=0.5*(1.0+1.0/eps);
// 	for (i=0; i<nface; i++){
// 		z[i]=r[i]/scale1;
// 		z[i+nface]=r[i+nface]/scale2;
// 	}
// 	return 0;

// }




/**********************************************************/
/* lapack provide lu decomposition, however, something    */
/* is wrong with cmake ************************************/
/**********************************************************/
int lu_decomp( double **A, int N, int *ipiv ) {

	int i, j, k, imax;
	double maxA, *ptr, absA, Tol = 1.0e-14;

  	for ( i = 0; i <= N; i++ )
   	ipiv[i] = i; // record pivoting number

  	for ( i = 0; i < N; i++ ) {
   	maxA = 0.0;
   	imax = i;
   	for (k = i; k < N; k++)
   	  	if ((absA = fabs(A[k][i])) > maxA) {
   	   	maxA = absA;
   	    	imax = k;
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
   	  	for (k = i + 1; k < N; k++)
   	   A[j][k] -= A[j][i] * A[i][k];
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

   	for (k = 0; k < i; k++)
      	xtemp[i] -= matrixA[i][k] * xtemp[k];
  	}

  	for (i = N - 1; i >= 0; i--) {
    	for (k = i + 1; k < N; k++)
      	xtemp[i] -= matrixA[i][k] * xtemp[k];

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
  
  	clock_t start_p,finish_p;
  	double total_p = 0;
  	start_p = clock();
	
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
    	leaflength(troot, idx);
    	nrow  = Nrow;
    	nrow2 = nrow*2;
    	ibeg  = idx;
    	iend  = idx + nrow - 1;

    	for ( i = ibeg; i <= iend; i++ ) {
      	tp[0] = tr_xyz[3*i];//s_particle_position[0][i];
      	tp[1] = tr_xyz[3*i+1]; //s_particle_position[1][i];
      	tp[2] = tr_xyz[3*i+2];//s_particle_position[2][i];
      	tq[0] = tr_q[3*i]; //s_particle_normal[0][i];
      	tq[1] = tr_q[3*i+1]; //s_particle_normal[1][i];
      	tq[2] = tr_q[3*i+2];//s_particle_normal[2][i];


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
  	finish_p = clock();
  	total_p = (double)(finish_p - start_p);
  	printf("psolve time is %f\n", total_p);
  	// return 0;

}
