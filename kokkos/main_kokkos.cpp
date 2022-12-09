/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
/* c */
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
#include "pp_timer.h"

/* c++ */
#include <cstdlib>
#include <cstdio>
#include <cmath>

/* kokkos */
// #include <sys/time.h>
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
// extern double *dev_xp, *dev_yp, *dev_zp, *dev_q, *dev_pot;

const double eps = 80.0;

#ifdef __cplusplus
extern "C"
{
#endif
// int *matvec ();
// int *psolve ();
// double eps;
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
	double s[3], pot=0.0, sum=0.0, pot_temp=0.0;	// potential
	double ptl, soleng, t1, t2;
	char fname[16], density[16];
	void readin(char fname[16], char density[16]);

	/*GMRES related variables*/
	static long int info;
	long int RESTRT, ldw, ldh, iter, N;
	double resid;
	extern void comp_source_wrapper();				// yang
	extern void comp_soleng_wrapper(double soleng);	// yang
	extern int *matvec(double *alpha, double *x, double *beta, double *y); // yang
	extern int *psolve(double *z, double *r); // yang
	// extern int *matvec(),*psolve();
	// extern int gmres_(long int *n, double *b, double *x, long int *restrt, double *work, long int *ldw, 
	// 	double *h, long int *ldh, long int *iter, double *resid, int *matvec (), int *psolve (), long int *info);

   extern void timer_start(char *n); // yang
   extern void timer_end(void); // yang

   Kokkos::initialize(argc, argv);
   {

   typedef Kokkos::Serial   HostExecSpace;
   typedef Kokkos::Cuda     DevExecSpace;
   typedef Kokkos::CudaSpace    MemSpace;
   typedef Kokkos::LayoutRight  Layout;
   typedef Kokkos::RangePolicy<HostExecSpace>  host_range_policy;
  	typedef Kokkos::RangePolicy<DevExecSpace>   dev_range_policy;


	timer_start((char*) "TOTAL_TIME");
	printf("%d %s %s %s \n", argc, argv[0], argv[1], argv[2]);

	/* read in structural information */
   sprintf(fname, "1a63");
   // sprintf(density, "1");
   // sprintf(fname,"%s",argv[1]);
   sprintf(density,"%s",argv[1]);
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

	// Allocate y, x vectors and Matrix A on device.
   // typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorType;
   // typedef Kokkos::View<double**, Layout, MemSpace>  ViewMatrixType;
   // ViewVectorType y( "y", N );
   // ViewVectorType x( "x", M );
   // ViewMatrixType A( "A", N, M );

	// xvct=(double *) calloc(N, sizeof(double));
	// work=(double *) calloc (ldw*(RESTRT+4), sizeof(double));
	// h=(double *) calloc (ldh*(RESTRT+2), sizeof(double));

	xvct=(double *) (Kokkos::kokkos_malloc(N * sizeof(double)));
	work=(double *) (Kokkos::kokkos_malloc(ldw*(RESTRT+4) * sizeof(double)));
	h=(double *) (Kokkos::kokkos_malloc(ldh*(RESTRT+2) * sizeof(double)));

   // Create host mirrors of device views.
   // ViewVectorType::HostMirror h_y = Kokkos::create_mirror_view( y );
   // ViewVectorType::HostMirror h_x = Kokkos::create_mirror_view( x );
   // ViewMatrixType::HostMirror h_A = Kokkos::create_mirror_view( A );

	gmres_(&N, bvct, xvct, &RESTRT, work, &ldw, h, &ldh, &iter, &resid, &matvec, &psolve, &info);

	Kokkos::fence();
	soleng=0.0;

	comp_soleng_wrapper(soleng); //wraps the solvation energy computation
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

	// free(tr_xyz);
	// free(tr_q);

	// free(tr_area);
	// free(bvct);
	// free(xvct);
	// free(atmchr);
	// free(chrpos);

	Kokkos::kokkos_free(tr_xyz);
	Kokkos::kokkos_free(tr_q);

	Kokkos::kokkos_free(tr_area);
	Kokkos::kokkos_free(bvct);
	Kokkos::kokkos_free(xvct);


  	Kokkos::kokkos_free(atmchr);
  	Kokkos::kokkos_free(chrpos);

	}
	Kokkos::finalize();


   return 0;
}

// ****************************************************************
int *psolve(double *z, double *r) {
/*r as original while z as scaled*/
	int i;
	double scale1,scale2;
	scale1=0.5*(1.0+eps);
	scale2=0.5*(1.0+1.0/eps);
	for (i=0; i<nface; i++){
		z[i]=r[i]/scale1;
		z[i+nface]=r[i+nface]/scale2;
	}
	return 0;

}
