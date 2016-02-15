/* BOUNDARY CONDITIONS
 *  Linear boundary conditions for the ODEs:
 *   r = A*s0 + B*s1 - c = 0
 *
 * REFERENCES
 *  Stoer & Bulirsch, pp. 557-561
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  par		[input]		parameters
 *  A			[output]	matrix of coefficients at t = 0
 *  B			[output]	matrix of coefficients at t = 1
 *  c			[output]	right-hand side vector
 */

#ifndef BCOND_H
#define BCOND_H


/* HEADER FILES */
#include <math.h>

/* PROTOTYPES */
void bcond(int, double*, double*, double*, double*);

/* IMPLEMENTATIONS */
/* Calculate coefficients for linear BC
 *   r = A*y(0) + B*y(1) - c = 0 */
void bcond(int n, double *par,
           double *A, double *B, double *c){
  int i; // for (i+1)st BC
  int j; // for jth element of y
  double PIH = M_PI/2.0;
  
  // initialize
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      A[i*n + j] = 0;
      B[i*n + j] = 0;
    }
    c[i] = 0;
  }
	
	// define parameters
	double v    = par[0];
	double area = 4.0*M_PI;
	double vlme = 4.0*M_PI*v/3.0;
  
	// matrix of coefficients for BC at t = 0
  A[0 *n + 0 ] =  1; // r  (0) = 0
  A[1 *n + 1 ] =  1; // psi(0) = -M_PI/2
  A[2 *n + 2 ] =  1; // p  (0) - p(1) = 1
  A[3 *n + 4 ] =  1; // A  (0) = 0
  A[4 *n + 5 ] =  1; // V  (0) = 0
  A[5 *n + 7 ] =  1; // R  (0) = 1
  A[6 *n + 11] =  1; // xcm(0) = 0
  
  // matrix of coefficients for BC at t = 1
  B[2 *n + 2 ] = -1; // p  (0) - p(1) = 1
  B[7 *n + 0 ] =  1; // r  (1) = 0
  B[8 *n + 1 ] =  1; // psi(1) = M_PI/2
  B[9 *n + 4 ] =  1; // A  (1) = area
  B[10*n + 5 ] =  1; // V  (1) = vlme
  B[11*n + 11] =  1; // xcm(1) = 0

  // right-hand side vector
	double a = sqrt(area/(4.0*M_PI));
  c[1 ] = -PIH;
  c[2 ] =  1.0;
  c[5 ] =  1.0/a;
  c[8 ] =  PIH;
  c[9 ] =  area;
  c[10] =  vlme;
}


#endif
