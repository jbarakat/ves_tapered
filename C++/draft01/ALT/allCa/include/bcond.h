/* BOUNDARY CONDITIONS
 *  Linear boundary conditions for the ODEs:
 *   r = A*s0 + B*s1 - c = 0
 *
 * REFERENCES
 *  Stoer & Bulirsch, pp. 557-561
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  par   [input]   parameters
 *  A			[output]	matrix of coefficients at t = 0
 *  B			[output]	matrix of coefficients at t = 1
 *  c			[output]	right-hand side vector
 */

#ifndef BCOND_H
#define BCOND_H


/* HEADER FILES */
#include <math.h>
#include <gsl/gsl_sf_trig.h>

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
	
	// parameters
	double T0   = par[0]; // membrane relaxation time = mu*a^3/kappa
	double alph = par[1]; // taper angle of tube wall
	double v    = par[2]; // reduced volume
	double dp   = par[3]; // pressure drop (scaled by kappa/a^3)
	double R0   = par[4]; // tube radius at center-of-mass axial position (scaled by a)
	double xcom = par[5]; // center-of-mass position 
	                      //   = time integral of center-of-mass translational speed

	double tana = gsl_sf_sin(alph)/gsl_sf_cos(alph);
  
  // initialize
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      A[i*n + j] = 0;
      B[i*n + j] = 0;
    }
    c[i] = 0;
  }

  // matrix of coefficients for BC at t = 0
  A[0 *n + 0 ] =  1   ; // r  (0) = 0
  A[1 *n + 2 ] =  1   ; // psi(0) = -PIH
  A[2 *n + 4 ] =  1   ; // qs (0) = 0
  A[3 *n + 7 ] =  1   ; // A  (0) = 0
  A[4 *n + 8 ] =  1   ; // V  (0) = 0
  A[5 *n + 10] =  1   ; // R  (0) + x(0)*tana = R0
  A[5 *n + 1 ] =  tana; 
  A[6 *n + 11] =  1   ; // X  (0) = 0
  A[7 *n + 5 ] =  1   ; // p  (0) - p(1) = dp
  
  // matrix of coefficients for BC at t = 1
  B[7 *n + 5 ] = -1; // p  (0) - p(1) = dp
  B[8 *n + 0 ] =  1; // r  (1) = 0
  B[9 *n + 2 ] =  1; // psi(1) = PIH
  B[10*n + 4 ] =  1; // qs (1) = 0 
  B[11*n + 7 ] =  1; // A  (1) = 4.0*M_PI
  B[12*n + 8 ] =  1; // V  (1) = 4.0*M_PI*v/3.0
  B[13*n + 11] =  1; // X  (1) = xcom

  // right-hand side vector
  c[2 ] = -PIH           ;
  c[5 ] =  R0            ;
  c[9 ] =  PIH           ;
  c[11] =  4.0*M_PI      ;
  c[12] =  4.0*M_PI*v/3.0;
}


#endif
