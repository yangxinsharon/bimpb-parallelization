/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
/* c */
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>
#include "gl_variables.h"
#include "gl_constants.h"


/* c++ */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <Kokkos_Core.hpp>


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
	double pre1, pre2;
    pre1=0.50*(1.0+eps); /* const eps=80.0 */
    pre2=0.50*(1.0+1.0/eps);
    Kokkos::parallel_for("matvecmul", nface, KOKKOS_LAMBDA(int i) {
    	double tp[3] = {tr_xyz[3*i], tr_xyz[3*i+1], tr_xyz[3*i+2]};
		double tq[3] = {tr_q[3*i], tr_q[3*i+1], tr_q[3*i+2]};

		double peng[2] = {0.0, 0.0};
		int j;
		for (j=0; j<nface; j++) {
        	if (j != i) {
				double sp[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
				double sq[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
				double r_s[3] = {sp[0]-tp[0], sp[1]-tp[1], sp[2]-tp[2]};
				double sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
				double rs = sqrt(sumrs);
				double irs = 1.0/sqrt(sumrs) ; //rsqrt(sumrs);
				double G0 = one_over_4pi;
				G0 = G0*irs;
				double kappa_rs = kappa*rs;
				double exp_kappa_rs = exp(-kappa_rs);
				double Gk = exp_kappa_rs*G0;
	
				double cos_theta = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
				double cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
	
				double tp1 = G0*irs;
				double tp2 = (1.0+kappa_rs) * exp_kappa_rs;
	
				double G10 = cos_theta0*tp1;
				double G20 = tp2*G10;
	
				double G1 = cos_theta*tp1;
				double G2 = tp2*G1;
	
				double dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
				double G3 = (dot_tqsq-3.0*cos_theta0*cos_theta) * irs*tp1;
				double G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
				double L1 = G1-eps*G2;							// K2
				double L2 = G0-Gk;								// K1
				double L3 = G4-G3;								// K4
				double L4 = G10-G20/eps; //fdivide(G20,eps);	// K3
	
				double peng_old[2] = {x[j], x[j+nface]};
				double area = tr_area[j];
				peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
				peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
        	}
		}

		y[i] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
		y[nface+i] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
	});

	// int i, j;
	// double pre1, pre2;
	// double area, rs, irs, sumrs;
	// double G0, kappa_rs, exp_kappa_rs, Gk;
	// double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
	// double G10, G20, G1, G2, G3, G4;
	// double L1, L2, L3, L4;

    // pre1=0.50*(1.0+eps); /* const eps=80.0 */
    // pre2=0.50*(1.0+1.0/eps);
    // for (i=0; i<nface; i++) {
    // 	double tp[3] = {tr_xyz[3*i], tr_xyz[3*i+1], tr_xyz[3*i+2]};
	// 	double tq[3] = {tr_q[3*i], tr_q[3*i+1], tr_q[3*i+2]};

	// 	double peng[2] = {0.0, 0.0};
	// 	for (j=0; j<nface; j++) {
    //     	if (j != i) {
	// 			double sp[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
	// 			double sq[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
	// 			double r_s[3] = {sp[0]-tp[0], sp[1]-tp[1], sp[2]-tp[2]};
	// 			sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
	// 			rs = sqrt(sumrs);
	// 			irs = 1.0/sqrt(sumrs) ; //rsqrt(sumrs);
	// 			G0 = one_over_4pi;
	// 			G0 = G0*irs;
	// 			kappa_rs = kappa*rs;
	// 			exp_kappa_rs = exp(-kappa_rs);
	// 			Gk = exp_kappa_rs*G0;
	
	// 			cos_theta = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
	// 			cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
	
	// 			tp1 = G0*irs;
	// 			tp2 = (1.0+kappa_rs) * exp_kappa_rs;
	
	// 			G10 = cos_theta0*tp1;
	// 			G20 = tp2*G10;
	
	// 			G1 = cos_theta*tp1;
	// 			G2 = tp2*G1;
	
	// 			dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
	// 			G3 = (dot_tqsq-3.0*cos_theta0*cos_theta) * irs*tp1;
	// 			G4 = tp2*G3 - kappa2*cos_theta0*cos_theta*Gk;
	// 			L1 = G1-eps*G2;							// K2
	// 			L2 = G0-Gk;								// K1
	// 			L3 = G4-G3;								// K4
	// 			L4 = G10-G20/eps; //fdivide(G20,eps);	// K3
	
	// 			double peng_old[2] = {x[j], x[j+nface]};
	// 			area = tr_area[j];
	// 			peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
	// 			peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
    //     	}
	// 	}

	// 	y[i] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
	// 	y[nface+i] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
	// }
    Kokkos::fence();
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

	// if ((chrptl=(double *) malloc(nface*sizeof(double)))==NULL) {
    if ((chrptl=(double *) (Kokkos::kokkos_malloc(nface*sizeof(double))))==NULL) {
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
	// int i, j;
    // double sumrs, irs, rs, G0, Gk, kappa_rs, exp_kappa_rs;
    // double cos_theta, G1, G2, L1, L2, tp1, tp2;
	// for (j=0; j<nface; j++) {
    // 	ptl[j] = 0.0;
	// 	double r[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
	// 	double v[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
    // 	for (i=0; i<nchr; i++) {
    //     	double s[3] = {chrpos[3*i], chrpos[3*i+1], chrpos[3*i+2]};
	// 		double r_s[3] = {r[0]-s[0], r[1]-s[1], r[2]-s[2]};
	// 		sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
	// 		rs = sqrt(sumrs);
	// 		irs = 1.0/sqrt(sumrs);

    //     	G0 = one_over_4pi;
    //     	G0 = G0*irs;
    //     	kappa_rs = kappa*rs;
    //     	exp_kappa_rs = exp(-kappa_rs);
    //     	Gk = exp_kappa_rs*G0;

    //     	cos_theta = (v[0]*r_s[0]+v[1]*r_s[1]+v[2]*r_s[2]) * irs;

    //     	tp1 = G0*irs;
    //     	tp2 = (1.0+kappa_rs)*exp_kappa_rs;

    //     	G1 = cos_theta*tp1;
    //     	G2 = tp2*G1;

    //     	L1 = G1-eps*G2;
    //     	L2 = G0-Gk;

    //   		ptl[j] = ptl[j] + atmchr[i] * (L1*xvct[j]+L2*xvct[nface+j]) * tr_area[j];
	// 	}
    // }
    Kokkos::parallel_for("comp_pot", nface, KOKKOS_LAMBDA(int j) {
    	ptl[j] = 0.0;
		double r[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
		double v[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
    	for (int i=0; i<nchr; i++) {
    	   	double s[3] = {chrpos[3*i], chrpos[3*i+1], chrpos[3*i+2]};
			double r_s[3] = {r[0]-s[0], r[1]-s[1], r[2]-s[2]};
			double sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
			double rs = sqrt(sumrs);
			double irs = 1.0/sqrt(sumrs);

    	   	double G0 = one_over_4pi;
    	   	G0 = G0*irs;
    	   	double kappa_rs = kappa*rs;
    	   	double exp_kappa_rs = exp(-kappa_rs);
    	   	double Gk = exp_kappa_rs*G0;

    	   	double cos_theta = (v[0]*r_s[0]+v[1]*r_s[1]+v[2]*r_s[2]) * irs;

    	   	double tp1 = G0*irs;
    	   	double tp2 = (1.0+kappa_rs)*exp_kappa_rs;

    	   	double G1 = cos_theta*tp1;
    	   	double G2 = tp2*G1;

    	   	double L1 = G1-eps*G2;
    	   	double L2 = G0-Gk;

    	  	ptl[j] = ptl[j] + atmchr[i] * (L1*xvct[j]+L2*xvct[nface+j]) * tr_area[j];
		}
    });
    Kokkos::fence();
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
	// int i, j;
	// double sumrs, cos_theta, irs, G0, G1, tp1;
	// for (i=0; i<nface; i++) {
    //     bvct[i] = 0.0;
    //     bvct[i+nface] = 0.0;
    //     for (j=0; j<nchr; j++) {
    //         double r_s[3] = {chrpos[3*j]-tr_xyz[3*i], chrpos[3*j+1]-tr_xyz[3*i+1], 
    //         	chrpos[3*j+2]-tr_xyz[3*i+2]};
	// 		sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2]; 
    //         cos_theta = tr_q[3*i]*r_s[0] + tr_q[3*i+1]*r_s[1] + tr_q[3*i+2]*r_s[2];
	// 		irs = 1.0/sqrt(sumrs) ;//rsqrt(sumrs);//returns reciprocal square root of scalars and vectors.
    //         cos_theta = cos_theta*irs;
    //         G0 = one_over_4pi;//constant
    //         G0 = G0*irs;
    //         tp1 = G0*irs;
    //         G1 = cos_theta*tp1;
    //         bvct[i] = bvct[i]+atmchr[j]*G0;
    //         bvct[nface+i] = bvct[nface+i]+atmchr[j]*G1;
    //     }

    // }
	Kokkos::parallel_for("comp_source", nface, KOKKOS_LAMBDA(int i) {
    	bvct[i] = 0.0;
    	bvct[i+nface] = 0.0;
    	for (int j=0; j<nchr; j++) {
    	    double r_s[3] = {chrpos[3*j]-tr_xyz[3*i], chrpos[3*j+1]-tr_xyz[3*i+1], 
    	    	chrpos[3*j+2]-tr_xyz[3*i+2]};
			double sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2]; 
    	    double cos_theta = tr_q[3*i]*r_s[0] + tr_q[3*i+1]*r_s[1] + tr_q[3*i+2]*r_s[2];
			double irs = 1.0/sqrt(sumrs) ;//rsqrt(sumrs);//returns reciprocal square root of scalars and vectors.
    	    cos_theta = cos_theta*irs;
    	    double G0 = one_over_4pi;//constant
    	    G0 = G0*irs;
    	    double tp1 = G0*irs;
    	    double G1 = cos_theta*tp1;
    	    bvct[i] = bvct[i]+atmchr[j]*G0;
    	    bvct[nface+i] = bvct[nface+i]+atmchr[j]*G1;
    	}
    });
    Kokkos::fence();
}
