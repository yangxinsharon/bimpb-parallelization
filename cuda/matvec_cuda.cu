#include <stdio.h>
#include "gl_variables.h"
#include "gl_constants.h"

extern "C" int *matvec(double *alpha, double *x, double *beta, double *y);
extern "C" void comp_soleng_wrapper(double soleng);
extern "C" void comp_source_wrapper();
extern "C" void initGPU();
extern "C" void freeGPU();
__global__ void comp_pot(const double* xvct, double *atmchr, double *chrpos,
double *ptl, double *tr_xyz,double *tr_q, double *tr_area, int nface, int nchr);
__global__ void comp_source( double* bvct, double *atmchr, double *chrpos,
double *tr_xyz, double *tr_q, int nface, int nchr);

#define checkcudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
    if(cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
        exit(-1);
    }
}
// This will output the proper error string when calling cudaGetLastError
#define getLastCudaError(msg)      __getLastCudaError (msg, __FILE__, __LINE__)

inline void __getLastCudaError(const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if (cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
        file, line, errorMessage, (int)err, cudaGetErrorString( err ) );
        exit(-1);
    }
}

// get the ID from block
__global__ void initY( double* y, int nface) {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        if(i<nface)
         y[i]=0.0;
}


__global__ void matvecmul(const double *x, double *y, double *q, int nface
        ,double *tr_xyz,double *tr_q, double *tr_area,double alpha, double beta){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j;
	double pre1,pre2;
	double area, rs, irs,sumrs;
	double G0, kappa_rs, exp_kappa_rs, Gk;
	double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
	double G10, G20, G1, G2, G3, G4;
	double L1, L2, L3, L4;

    pre1=0.50*(1.0+eps); /* eps=80.0 a constant */
    pre2=0.50*(1.0+fdivide(1.0,eps)); //!!! fdivide?
    if (i < nface) {

		double3 tp={tr_xyz[3*i],tr_xyz[3*i+1],tr_xyz[3*i+2]};
		double3 tq={tr_q[3*i],tr_q[3*i+1],tr_q[3*i+2]};

		double2 peng={0.0,0.0};
		for (j=0;j<nface;j++){
          if (j != i){
			double3 sp={tr_xyz[3*j],tr_xyz[3*j+1],tr_xyz[3*j+2]};
			double3 sq={tr_q[3*j],tr_q[3*j+1],tr_q[3*j+2]};
			double3 r_s={sp.x-tp.x,sp.y-tp.y,sp.z-tp.z};
			sumrs= r_s.x*r_s.x + r_s.y*r_s.y+r_s.z*r_s.z;
			rs=sqrt(sumrs);
			irs=rsqrt(sumrs);
			G0=one_over_4pi;
			G0=G0*irs;
			kappa_rs=kappa*rs;
			exp_kappa_rs=exp(-kappa_rs);
			Gk=exp_kappa_rs*G0;

			cos_theta	=(sq.x*r_s.x+sq.y*r_s.y+sq.z*r_s.z)*irs;
			cos_theta0	=(tq.x*r_s.x+tq.y*r_s.y+tq.z*r_s.z)*irs;

			tp1=G0*irs;
			tp2=(1.0+kappa_rs)*exp_kappa_rs;

			G10=cos_theta0*tp1;
			G20=tp2*G10;

			G1=cos_theta*tp1;
			G2=tp2*G1;

			dot_tqsq=sq.x*tq.x+sq.y*tq.y+sq.z*tq.z;
			G3=(dot_tqsq-3.0*cos_theta0*cos_theta)*irs*tp1;
			G4=tp2*G3-kappa2*cos_theta0*cos_theta*Gk;
			L1=G1-eps*G2;
			L2=G0-Gk;
			L3=G4-G3;
			L4=G10-fdivide(G20,eps);

			double2 peng_old={x[j], x[j+nface]};
			area=tr_area[j];
			peng.x=peng.x+(L1*peng_old.x+L2*peng_old.y)*area;
			peng.y=peng.y+(L3*peng_old.x+L4*peng_old.y)*area;
          }
		}

		y[i]=y[i] *beta+(pre1*x[i]-peng.x) * alpha;
		y[nface+i]=y[nface+i] * beta+(pre2*x[nface+i]-peng.y) * alpha;
	}
}

double *d_X, *d_Y,*d_tr_xyz,*d_tr_q,*d_tr_area,*d_atmchr,*d_chrpos,*d_chrptl,*d_xvct;
int threadsPerBlock = 256;

void initGPU() {
	checkcudaErrors(cudaMalloc((void**)&d_X,2*nface*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_Y,2*nface*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_tr_q,3*nface*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_tr_xyz,3*nface*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_tr_area,nface*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_atmchr,nchr*sizeof(double))) ;
    checkcudaErrors(cudaMalloc((void**)&d_chrpos,3*nchr*sizeof(double))) ;

    checkcudaErrors(cudaMemcpy(d_tr_area,tr_area,nface*sizeof(double),cudaMemcpyHostToDevice));
    checkcudaErrors(cudaMemcpy(d_tr_xyz,tr_xyz,3*nface*sizeof(double),cudaMemcpyHostToDevice));
    checkcudaErrors(cudaMemcpy(d_tr_q,tr_q,3*nface*sizeof(double),cudaMemcpyHostToDevice));
    checkcudaErrors(cudaMemcpy(d_atmchr,atmchr,nchr*sizeof(double),cudaMemcpyHostToDevice));
    checkcudaErrors(cudaMemcpy(d_chrpos,chrpos,3*nchr*sizeof(double),cudaMemcpyHostToDevice));
}


void freeGPU() {
    checkcudaErrors( cudaFree(d_tr_area));
    checkcudaErrors( cudaFree(d_tr_xyz));
    checkcudaErrors( cudaFree(d_tr_q));
    checkcudaErrors( cudaFree(d_X));
    checkcudaErrors( cudaFree(d_Y));
    checkcudaErrors( cudaFree(d_xvct));
    checkcudaErrors( cudaFree(d_atmchr));
    checkcudaErrors( cudaFree(d_chrpos));
    checkcudaErrors( cudaFree(d_chrptl));
}

/* This subroutine wraps the matrix-vector multiplication */
int *matvec(double *alpha, double *x, double *beta, double *y) {
	int blocksPerGrid = (nface + threadsPerBlock - 1) / threadsPerBlock;
    checkcudaErrors(cudaMemcpy(d_X, x,2*nface*sizeof(double), cudaMemcpyHostToDevice));
    checkcudaErrors(cudaMemcpy(d_Y, y,2*nface*sizeof(double), cudaMemcpyHostToDevice));
    matvecmul<<<blocksPerGrid, threadsPerBlock>>>(d_X, d_Y,d_tr_q,nface
            ,d_tr_xyz,d_tr_q,d_tr_area, *alpha, *beta);
    getLastCudaError("kernel launch failure");
    checkcudaErrors(cudaMemcpy(y, d_Y, 2*nface*sizeof(double), cudaMemcpyDeviceToHost));
    return NULL;
}

/* This subroutine wraps the solvation energy computation */
/* Called before freeGPU() */
void comp_soleng_wrapper(double soleng) {
    int i;
	double *chrptl;
	int blocksPerGrid = (nface + threadsPerBlock - 1) / threadsPerBlock;
	double units_para=2.0;
    units_para=units_para*units_coef;
    units_para=units_para*pi;

	if ((chrptl=(double *) malloc(nface*sizeof(double)))==NULL){
		printf("error in allcating chrptl");
	}
	checkcudaErrors(cudaMalloc((void**)&d_chrptl,nface*sizeof(double))) ;
	checkcudaErrors(cudaMalloc((void**)&d_xvct,2*nface*sizeof(double))) ;
    checkcudaErrors(cudaMemcpy(d_xvct,xvct,2*nface*sizeof(double),cudaMemcpyHostToDevice));

	comp_pot<<<blocksPerGrid, threadsPerBlock>>>(d_xvct, d_atmchr, d_chrpos,
    d_chrptl, d_tr_xyz,d_tr_q, d_tr_area, nface, nchr);
    checkcudaErrors(cudaMemcpy(chrptl,d_chrptl,nface*sizeof(double),cudaMemcpyDeviceToHost));

	soleng=0.0;
	for (i=0;i<nface;i++) soleng=soleng+chrptl[i];
	soleng=soleng*units_para;
	printf("solvation energy on GPU = %f kcal/mol\n",soleng);
}



/* This subroutine calculates the element-wise potential on GPU */
__global__ void comp_pot(const double* xvct, double *atmchr, double *chrpos,
double *ptl, double *tr_xyz,double *tr_q, double *tr_area, int nface, int nchr){

	int j = blockDim.x * blockIdx.x + threadIdx.x;
    double sumrs,irs,rs,G0,Gk,kappa_rs,exp_kappa_rs;
    double cos_theta,G1,G2,L1,L2,tp1,tp2;
    int i;
	if (j<nface){
    	ptl[j]=0.0;
		double3 r={tr_xyz[3*j],tr_xyz[3*j+1],tr_xyz[3*j+2]};
		double3 v={tr_q[3*j],tr_q[3*j+1],tr_q[3*j+2]};
    	for (i=0; i<nchr; i++) {
        	double3 s={chrpos[3*i],chrpos[3*i+1],chrpos[3*i+2]};
			double3 r_s={r.x-s.x,r.y-s.y,r.z-s.z};
			sumrs= r_s.x*r_s.x + r_s.y*r_s.y+r_s.z*r_s.z;
			rs=sqrt(sumrs);
			irs=rsqrt(sumrs);

        	G0=one_over_4pi;
        	G0=G0*irs;
        	kappa_rs=kappa*rs;
        	exp_kappa_rs=exp(-kappa_rs);
        	Gk=exp_kappa_rs*G0;

        	cos_theta=(v.x*r_s.x+v.y*r_s.y+v.z*r_s.z)*irs;

        	tp1=G0*irs;
        	tp2=(1.0+kappa_rs)*exp_kappa_rs;

        	G1=cos_theta*tp1;
        	G2=tp2*G1;

        	L1=G1-eps*G2;
        	L2=G0-Gk;

      		ptl[j]=ptl[j]+atmchr[i]*(L1*xvct[j]+L2*xvct[nface+j])*tr_area[j];
		}
    }
}

/* This subroutine wraps the solvation energy computation */
/* In main_cuda.c after initGPU() */
void comp_source_wrapper() {
	double *d_bvct;
    int blocksPerGrid = (nface + threadsPerBlock - 1) / threadsPerBlock;

	checkcudaErrors(cudaMalloc((void**)&d_bvct,2*nface*sizeof(double))) ;

	comp_source<<<blocksPerGrid, threadsPerBlock>>>(d_bvct, d_atmchr, d_chrpos,
    d_tr_xyz,d_tr_q, nface, nchr);
    checkcudaErrors(cudaMemcpy(bvct,d_bvct,2*nface*sizeof(double),cudaMemcpyDeviceToHost));
    checkcudaErrors(cudaFree(d_bvct));
}


/* This subroutine calculates the source term of the integral equation on GPU */
/* atmchr=atom charge   chrpos=charge position */
/* bvct be located at readin.c */
__global__ void comp_source( double* bvct, double *atmchr, double *chrpos,
double *tr_xyz,double *tr_q, int nface, int nchr){

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j;
	double sumrs,cos_theta,irs,G0,G1,tp1;
	if (i<nface){
        bvct[i]=0.0;
        bvct[i+nface]=0.0;
        for (j=0; j<nchr; j++){
            double3 r_s={	chrpos[3*j]-tr_xyz[3*i],
							chrpos[3*j+1]-tr_xyz[3*i+1],
							chrpos[3*j+2]-tr_xyz[3*i+2]};
			sumrs= r_s.x*r_s.x + r_s.y*r_s.y+r_s.z*r_s.z; //c can't use that r_s.x
            cos_theta=tr_q[3*i]*r_s.x+tr_q[3*i+1]*r_s.y+tr_q[3*i+2]*r_s.z;
			irs=rsqrt(sumrs);//returns reciprocal square root of scalars and vectors.
            cos_theta=cos_theta*irs;
            G0=one_over_4pi;//constant
            G0=G0*irs;
            tp1=G0*irs;
            G1=cos_theta*tp1;
            bvct[i]=bvct[i]+atmchr[j]*G0;
            bvct[nface+i]=bvct[nface+i]+atmchr[j]*G1;
        }

    }
}
