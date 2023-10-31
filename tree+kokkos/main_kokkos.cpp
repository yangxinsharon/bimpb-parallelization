/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include "pp_timer.h"

/* c++ */
#include <cstdlib>
#include <cstdio>
#include <cmath>
// #include "gl_constants.h"

/* kokkos */
#include <Kokkos_Core.hpp>

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
// extern double *dev_tr_xyz, *dev_tr_q, *dev_tr_area, *dev_bvct;
extern const double eps;
extern double **tr_xyz2D, **tr_q2D;
// extern int *ipiv;
// extern double *rhs;
// extern double **matrixA;
extern int *leafarr;
extern double *xtemp;
extern int Nleaf;
// extern int arridx;
// extern TreeNode *s_tree_root;

#ifdef __cplusplus
extern "C"
{
#endif
int gmres_(long int *n, double *b, double *x, long int *restrt, double *work, long int *ldw, 
		double *h, long int *ldh, long int *iter, double *resid, 
		int *matvec (double *alpha, double *x, double *beta, double *y), 
		int *psolve (double *z, double *r), long int *info);

#ifdef __cplusplus
}
#endif

int main(int argc, char *argv[]) {
	/*variables local to main*/
	int i,j;
	// double s[3], pot=0.0, sum=0.0, pot_temp=0.0;	// potential
	double ptl, soleng, t1, t2;
	char fname[16], density[16];
	extern double** Make2DDoubleArray(int arraySizeX, int arraySizeY, char info[]);
	extern void readin(char fname[16], char density[16]);

	/*GMRES related variables*/
	static long int info;
	long int RESTRT, ldw, ldh, iter, N;
	double resid;
	extern void comp_source_wrapper();				// yang
	extern void comp_soleng_wrapper(double soleng);	// yang
	extern int *matvec(double *alpha, double *x, double *beta, double *y); // yang
	extern int *psolve(double *z, double *r); // yang

   extern void timer_start(char *n); // yang
   extern void timer_end(void); // yang
   extern int TreecodeInitialization();
	extern int TreecodeFinalization();
	// extern void leaflength(TreeNode *p, int idx);

   Kokkos::initialize(argc, argv);
   {

	timer_start((char*) "TOTAL_TIME");
	printf("%d %s %s %s \n", argc, argv[0], argv[1], argv[2]);

	/* read in structural information */
   // sprintf(fname, "1ajj");
   // sprintf(density, "1");
   sprintf(fname,"%s",argv[1]);
   sprintf(density,"%s",argv[2]);
	readin(fname, density);
	comp_source_wrapper(); //wraps the solvation energy computation
	Kokkos::fence();

	/* parameters for GMRES */
	RESTRT=10;
	N=2*nface;
	ldw=N;
	ldh=RESTRT+1;
	iter=100;
	resid=1e-4;

	xvct=(double *) (Kokkos::kokkos_malloc(N * sizeof(double)));
	work=(double *) (Kokkos::kokkos_malloc(ldw*(RESTRT+4) * sizeof(double)));
	h=(double *) (Kokkos::kokkos_malloc(ldh*(RESTRT+2) * sizeof(double)));


	TreecodeInitialization();
	// Kokkos::fence();

	// leafarr = (int *) Kokkos::kokkos_malloc(3*Nleaf* sizeof(int));
	// int idx = 0, nrow = 0, ibeg = 0, iend = 0;
	// arridx = 0; // extern variable
	// while ( idx < nface ) {
	//    leaflength(s_tree_root, idx);
	//    nrow  = Nrow;
	//    ibeg  = idx;
	//    iend  = idx + nrow - 1;	
	//    leafarr[0+3*arridx] = ibeg;
	//    leafarr[1+3*arridx] = nrow;
	//    leafarr[2+3*arridx] = iend;    
	// 	// printf("ibeg iend nrow is %d, %d, %d\n",ibeg,iend,nrow);
	// 	arridx += 1;
	// 	idx += nrow;
	// }

	gmres_(&N, bvct, xvct, &RESTRT, work, &ldw, h, &ldh, &iter, &resid, &matvec, &psolve, &info);
// int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
//            info)
	Kokkos::fence();
	soleng=0.0;

	comp_soleng_wrapper(soleng); //wraps the solvation energy computation
	Kokkos::fence();
	// timer_end();
	TreecodeFinalization();
	Kokkos::fence();
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

	for(i=0;i<3;i++) {
		free(tr_xyz2D[i]);
	}	
	free(tr_xyz2D);

	for(i=0;i<3;i++) {
		free(tr_q2D[i]);
	}
	free(tr_q2D);

	// for(i=0;i<2*maxparnode;i++) {
	// 	free(matrixA[i]);
	// }	
	// free(matrixA);
	// Kokkos::kokkos_free(leafarr);  

	Kokkos::kokkos_free(tr_xyz);
	Kokkos::kokkos_free(tr_q);

	Kokkos::kokkos_free(tr_area);
	Kokkos::kokkos_free(bvct);
	Kokkos::kokkos_free(xvct);


  	Kokkos::kokkos_free(atmchr);
  	Kokkos::kokkos_free(chrpos);


 	// *psolve
  	// Kokkos::kokkos_free(rhs);
	// Kokkos::kokkos_free(ipiv);
  	// Kokkos::kokkos_free(leafarr);
	}
	Kokkos::finalize();


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
