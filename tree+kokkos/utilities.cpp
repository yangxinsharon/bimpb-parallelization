/**************************************************************************
* FILE NAME: utilities.c                                                  *
*                                                                         *
* PURPOSE: Contains helper routines for calculating min and max values    *
*          in an array (used by treecode.c and tabipb.c), as well as      *
*          routine to calculate triangle area (used by readin.c and       *
*          tabipb.c)                                                      *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 01/15/2018  Leighton Wilson   Moved TriangleArea from readin.c, moved   *
*                               Partition to its own file                 *
* 01/12/2018  Leighton Wilson   Localized all global variables            *
*                                                                         *
**************************************************************************/

/* c */
// #include <math.h>

/* c++ */
#include <cmath>
#include "utilities.h"


/**********************************************************/
double MinVal(double *variables, int number) 
{
    int i;
    double min_val;

    min_val = variables[0];
    for (i = 1; i < number; i++) {
        if (min_val > variables[i]) {
            min_val = variables[i];
        }
    }
    
    return min_val;
}
/**********************************************************/


/**********************************************************/
double MaxVal(double *variables, int number) 
{
    int i;
    double max_val;
  
    max_val = variables[0];
    for (i = 1; i < number; i++) {
        if (max_val < variables[i])
            max_val = variables[i];
    }
    
    return max_val;
}
/**********************************************************/


/**********************************************************/
/* function computing the area of a triangle given vertices coodinates */
double TriangleArea(double v[3][3])
{
    int i;
    double a[3], b[3], c[3], aa, bb, cc, ss, area;

    for (i = 0; i <= 2; i++) {
        a[i] = v[i][0] - v[i][1];
        b[i] = v[i][0] - v[i][2];
        c[i] = v[i][1] - v[i][2];
    }

    aa = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    bb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    cc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    ss = 0.5 * (aa + bb + cc);
    area = sqrt(ss * (ss-aa) * (ss-bb) * (ss-cc));

    return area;
}
/**********************************************************/


/**********************************************************/
/* lapack provide lu decomposition, however, something    */
/* is wrong with cmake */
// KOKKOS_FUNCTION
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
/**********************************************************/


/**********************************************************/
// KOKKOS_FUNCTION
void lu_solve( double **matrixA, int N, int *ipiv, double *rhs ) {
    /* b will contain the solution */
    double *xtemp;

    // make_vector(xtemp, N);
    xtemp=(double *) calloc(N, sizeof(double));
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
    // free_vector(xtemp);
    free(xtemp);
}
/**********************************************************/


