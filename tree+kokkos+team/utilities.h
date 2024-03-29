/**************************************************************************
* FILE NAME: utilities.h                                                  *
*                                                                         *
* PURPOSE: header for helper routines                                     *
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
* 01/15/2018  Leighton Wilson   Moved TriangleArea from readin.c          *
* 01/12/2018  Leighton Wilson   Localized all global variables            *
*                                                                         *
**************************************************************************/

#ifndef H_UTILITY_ROUTINES_H
#define H_UTILITY_ROUTINES_H

#include <Kokkos_Core.hpp>
	
double MinVal(double *variables, int number);

double MaxVal(double *variables, int number);

double TriangleArea(double v[3][3]);

// KOKKOS_FUNCTION
int lu_decomp( double **A, int N, int *ipiv );

// KOKKOS_FUNCTION
void lu_solve( double **matrixA, int N, int *ipiv, double *rhs );

#endif /* H_UTILITY_ROUTINES_H */
