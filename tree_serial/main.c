/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pp_timer.h"

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

extern double **tr_xyz2D, **tr_q2D;

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
   extern int TreecodeInitialization();
	extern int TreecodeFinalization();
	
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
	// iter=100;
	iter=10000;
	resid=1e-4;
	xvct=(double *) calloc(N, sizeof(double));

	work=(double *) calloc (ldw*(RESTRT+4), sizeof(double));
	h=(double *) calloc (ldh*(RESTRT+2), sizeof(double));


	TreecodeInitialization();

	gmres_(&N, bvct, xvct, &RESTRT, work, &ldw, h, &ldh, &iter, &resid, &matvec, &psolve, &info);

	soleng=0.0;

	comp_soleng_wrapper(soleng); //wraps the solvation energy computation
	timer_end();


	TreecodeFinalization();

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

	for(i=0;i<3;i++) {
		free(tr_xyz2D[i]);
	}	
	free(tr_xyz2D);

	for(i=0;i<3;i++) {
		free(tr_q2D[i]);
	}
	free(tr_q2D);

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
