/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gl_variables.h"
#include "gl_constants.h"
#include "mpi.h"
#include <string.h>

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
	/* declarations for mpi */
	int is, ie, i, j;
	int ierr, numprocs, myid;

	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_size = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_rank = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

    double pre1=0.50*(1.0+eps); /* const eps=80.0 */
    double pre2=0.50*(1.0+1.0/eps);

	/* determine this processor's interval */
  	is = ((int) (1.0*nface/numprocs))*myid;
  	ie = ((int) (1.0*nface/numprocs))*(myid+1);
  	if (myid == numprocs-1)  ie = nface;

  	int chunk=ie-is;
  	int irecv[numprocs];
	for (i=0; i<numprocs-1; i++) {
		irecv[i] = ((int) (1.0*nface/numprocs));
	}
	irecv[numprocs-1]=nface-((int) (1.0*nface/numprocs))*(numprocs-1);


	int idisp[numprocs];
	int sumtmp=0;
	for (i = 0; i <numprocs; i++) {
		idisp[i] = 0 + sumtmp;
		sumtmp = sumtmp+irecv[i];
	}

  	double *sbuf_y1, *sbuf_y2;
  	sbuf_y1 = (double *) calloc(chunk, sizeof(double));
	sbuf_y2 = (double *) calloc(chunk, sizeof(double));

    for (i=is; i<ie; i++) { 

    	double tp[3] = {tr_xyz[3*i], tr_xyz[3*i+1], tr_xyz[3*i+2]};
		double tq[3] = {tr_q[3*i], tr_q[3*i+1], tr_q[3*i+2]};

		double peng[2] = {0.0, 0.0};
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
				double L1 = G1-eps*G2;								// K2
				double L2 = G0-Gk;									// K1
				double L3 = G4-G3;									// K4
				double L4 = G10-G20/eps; //fdivide(G20,eps);		// K3
	
				double peng_old[2] = {x[j], x[j+nface]};
				double area = tr_area[j];
				peng[0] = peng[0] + (L1*peng_old[0] + L2*peng_old[1]) * area;
				peng[1] = peng[1] + (L3*peng_old[0] + L4*peng_old[1]) * area;
        	}
		}
		// y[i] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
  		// y[nface+i] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
  		sbuf_y1[i-idisp[myid]] = y[i]*beta + (pre1*x[i]-peng[0])*alpha;
  		sbuf_y2[i-idisp[myid]] = y[nface+i]*beta + (pre2*x[nface+i]-peng[1])*alpha;
	}


	double *rbuf_y1,*rbuf_y2;
	rbuf_y1 = (double *) calloc(nface, sizeof(double));
	rbuf_y2 = (double *) calloc(nface, sizeof(double));



	ierr = MPI_Allgatherv(sbuf_y1, chunk, MPI_DOUBLE, rbuf_y1, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather1 = %i\n",ierr);
  	}

	ierr = MPI_Allgatherv(sbuf_y2, chunk, MPI_DOUBLE, rbuf_y2, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather2 = %i\n",ierr);
  	}


  	for (i=0; i<nface; i++) {
  		y[i] = rbuf_y1[i];
  		y[nface+i] = rbuf_y2[i];
  	}

  	free(sbuf_y1);
  	free(sbuf_y2);
  	free(rbuf_y1);
  	free(rbuf_y2);
}



/* This subroutine wraps the matrix-vector multiplication */
int *matvec(double *alpha, double *x, double *beta, double *y) {
    matvecmul(x, y, tr_q, nface, tr_xyz, tr_q, tr_area, *alpha, *beta);
    // return NULL;
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
    int js, je, ierr, numprocs, myid;

	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_size = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_rank = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

  	js = ((int) (1.0*nface/numprocs))*myid;
  	je = ((int) (1.0*nface/numprocs))*(myid+1);
  	if (myid == numprocs-1)  je = nface;

	int chunk=je-js;
  	int irecv[numprocs];
	for (i=0; i<numprocs-1; i++) {
		irecv[i] = ((int) (1.0*nface/numprocs));
	}
	irecv[numprocs-1]=nface-((int) (1.0*nface/numprocs))*(numprocs-1);
	
	int idisp[numprocs];
	int sumtmp=0;
	for (i = 0; i <numprocs; i++) {
		idisp[i] = 0 + sumtmp;
		sumtmp = sumtmp+irecv[i];
	}

  	double *sbuf_ptl;
  	sbuf_ptl = (double *) calloc(chunk, sizeof(double));



	// for (j=0; j<nface; j++) {
	for (j=js; j<je; j++) {	
    	ptl[j] = 0.0;
		double r[3] = {tr_xyz[3*j], tr_xyz[3*j+1], tr_xyz[3*j+2]};
		double v[3] = {tr_q[3*j], tr_q[3*j+1], tr_q[3*j+2]};
    	for (i=0; i<nchr; i++) {
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
		sbuf_ptl[j-idisp[myid]]= ptl[j];
    }

    double *rbuf_ptl;
	rbuf_ptl = (double *) calloc(nface, sizeof(double));


	ierr = MPI_Allgatherv(sbuf_ptl, chunk, MPI_DOUBLE, rbuf_ptl, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather1 = %i\n",ierr);
  	}

  	for (i=0; i<nface; i++) {
  		ptl[i] = rbuf_ptl[i];
  	}

  	free(sbuf_ptl);
  	free(rbuf_ptl);

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
	/* declarations for mpi */
	int is, ie, i, j;
	int ierr, numprocs, myid;

	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_size = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (ierr != 0) {
	  printf(" error in MPI_Comm_rank = %i\n",ierr);
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

	/* determine this processor's interval */
  	is = ((int) (1.0*nface/numprocs))*myid;
  	ie = ((int) (1.0*nface/numprocs))*(myid+1);
  	if (myid == numprocs-1)  ie = nface;

	int chunk=ie-is;
  	int irecv[numprocs];
	for (i=0; i<numprocs-1; i++) {
		irecv[i] = ((int) (1.0*nface/numprocs));
	}
	irecv[numprocs-1]=nface-((int) (1.0*nface/numprocs))*(numprocs-1);
	
	int idisp[numprocs];
	int sumtmp=0;
	for (i = 0; i <numprocs; i++) {
		idisp[i] = 0 + sumtmp;
		sumtmp = sumtmp+irecv[i];
	}

  	double *sbuf_bvct1, *sbuf_bvct2;
  	sbuf_bvct1 = (double *) calloc(chunk, sizeof(double));
	sbuf_bvct2 = (double *) calloc(chunk, sizeof(double));

	for (i=is; i<ie; i++) {
        bvct[i] = 0.0;
        bvct[i+nface] = 0.0;
        for (j=0; j<nchr; j++) {
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
        sbuf_bvct1[i-idisp[myid]] = bvct[i];
		sbuf_bvct2[i-idisp[myid]] = bvct[nface+i];
    }
    double *rbuf_bvct1,*rbuf_bvct2;
	rbuf_bvct1 = (double *) calloc(nface, sizeof(double));
	rbuf_bvct2 = (double *) calloc(nface, sizeof(double));

	ierr = MPI_Allgatherv(sbuf_bvct1, chunk, MPI_DOUBLE, rbuf_bvct1, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather1 = %i\n",ierr);
  	}

  	ierr = MPI_Allgatherv(sbuf_bvct2, chunk, MPI_DOUBLE, rbuf_bvct2, irecv, idisp, MPI_DOUBLE, MPI_COMM_WORLD);
  	if (ierr != MPI_SUCCESS) {
  	   	printf("Error in MPI_Allgather2 = %i\n",ierr);
  	}
  	for (i=0; i<nface; i++) {
  		bvct[i] = rbuf_bvct1[i];
  		bvct[nface+i] = rbuf_bvct2[i];
  	}

  	free(sbuf_bvct1);
  	free(sbuf_bvct2);
  	free(rbuf_bvct1);
  	free(rbuf_bvct2);
}
