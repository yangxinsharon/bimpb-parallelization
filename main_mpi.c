/* Xin (Sharon) Yang
   SMU Mathematics Math 6370 
   Project bim-pb with Weihua Geng, Jiahui Chen */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pp_timer.h"
#include "mpi.h"

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

  	// initialize MPI
  	int ierr, numprocs, myid;
  	ierr = MPI_Init(&argc, &argv);
  	if (ierr != MPI_SUCCESS) {
  	   printf("Error in MPI_Init = %i\n",ierr);
  	   return 1;
  	}

  	ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	if (ierr != MPI_SUCCESS) {
  	   printf("Error in MPI_Comm_size = %i\n",ierr);
  	   return 1;
  	}

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  	if (ierr != MPI_SUCCESS) {
  	   printf("Error in MPI_Comm_rank = %i\n",ierr);
  	   return 1;
  	}

  	// start timer
  	double stime = MPI_Wtime();
  	// bimpb timer
	timer_start("TOTAL_TIME"); 

  	// root calling readin.c
  	if (myid == 0) {

		printf("%d %s %s %s \n", argc, argv[0], argv[1], argv[2]);

		/* read in structural information */
   	sprintf(fname, "1ajj");
   	// sprintf(density, "1");
   	// sprintf(fname,"%s",argv[1]);
   	sprintf(density,"%s",argv[1]);
		readin(fname, density);

	}
		

	// broadcast nface,nspt,natm,nchr from readin
	ierr = MPI_Bcast(&nface, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	if (ierr != 0) {
        printf("Error in MPI_Bcast nface = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
     }
	ierr = MPI_Bcast(&nspt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	if (ierr != 0) {
        printf("Error in MPI_Bcast nspt = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}
	ierr = MPI_Bcast(&natm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	if (ierr != 0) {
        printf("Error in MPI_Bcast natm = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}
	ierr = MPI_Bcast(&nchr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    	if (ierr != 0) {
        printf("Error in MPI_Bcast nchr = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}

	// allocate addresses that has been allocated in readin.c for other processes
	if (myid != 0){
		atmchr=(double *) malloc(nchr*sizeof(double));
		chrpos=(double *) malloc(3*nchr*sizeof(double));
		tr_xyz=(double *) calloc(3*nface, sizeof(double));
		tr_q=(double *) calloc(3*nface, sizeof(double));
		tr_area=(double *) calloc(nface, sizeof(double));
		bvct=(double *) calloc(2*nface, sizeof(double));
	}
	ierr = MPI_Barrier(MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Barrier = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}

	ierr = MPI_Bcast(atmchr, nchr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Bcast atmchr = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}

	ierr = MPI_Bcast(chrpos, 3*nchr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Bcast chrpos = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}		

	ierr = MPI_Bcast(tr_xyz, 3*nface, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Bcast tr_xyz = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}

	ierr = MPI_Bcast(tr_q, 3*nface, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Bcast tr_q = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}


	ierr = MPI_Bcast(tr_area, nface, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Bcast tr_area = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}

	ierr = MPI_Barrier(MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Barrier = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}	


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
	
	
	ierr = MPI_Barrier(MPI_COMM_WORLD);
     if (ierr != 0) {
        printf("Error in MPI_Barrier = %i\n",ierr);
        ierr = MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
	}
	
	timer_end();

	/* free memory */
	if (myid == 0){

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
	}
	else	{
		free(tr_xyz);
		free(tr_q);

		free(tr_area);
		free(bvct);
		free(xvct);
		free(atmchr);
		free(chrpos);
	}

	double ftime = MPI_Wtime();
	
	// finalize MPI
  	ierr = MPI_Finalize();

   return 0;
} // end main

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
