/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  Tozeren and Skalak, Int J Multiphase Flow (1979)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  par         [input]   parameters
 *  ur					[input]		rate of change of r
 *  t						[input]		abscissas
 *  y						[input]		solution
 *  f						[output]	function
 */

#ifndef FUNC_H
#define FUNC_H


/* HEADER FILES */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

/* PROTOTYPES */
void func(int, double*, double, double, double*, double*);

/* IMPLEMENTATIONS */
void func(int ny   , double *par, double ur,
          double t , double *y  , double *f){
	if (ny != 14)
    cout << "Error: ny should equal 14." << endl;

  int i;
  
	// characteristic scales:
	//  - a     = nominal radius of vesicle
	//  - kappa = bending modulus
	//  - mu    = suspending fluid viscosity

	// parameters
	double T0   = par[0]; // membrane relaxation time = mu*a^3/kappa
	double alph = par[1]; // taper angle of tube wall
	double v    = par[2]; // reduced volume
	double dp   = par[3]; // pressure drop (scaled by kappa/a^3)
	double R0   = par[4]; // tube radius at center-of-mass axial position (scaled by a)
	double xcom = par[5]; // center-of-mass position 
	                      //   = time integral of center-of-mass translational speed

	double tana = gsl_sf_sin(alph)/gsl_sf_cos(alph);
  
  // define variables from current timestep
  double r    = y[0 ]; // azimuthal radius of curvature		(scaled by a)
  double x    = y[1 ]; // axial position										(scaled by a)
  double psi  = y[2 ]; // tilt angle 
  double cs   = y[3 ]; // meridional curvature							(scaled by 1/a)
  double qs   = y[4 ]; // transverse shear tension         (scaled by kappa/a^2)
  double p    = y[5 ]; // gap pressure                     (scaled by kappa/a^3)
  double sig  = y[6 ]; // mean in-plane tension            (scaled by kappa/a^2)
  double A    = y[7 ]; // surface area                     (scaled by a^2)
  double V    = y[8 ]; // volume                           (scaled by a^3)
  double Q    = y[9 ]; // leakback flux per unit circumf.  (scaled by kappa/(mu*a))
  double R    = y[10]; // tube radius                      (scaled by a)
	double X    = y[11]; // center-of-mass axial position    (scaled by x)
  double U    = y[12]; // center-of-mass axial speed       (scaled by kappa/(mu*a^2))
  double S    = y[13]; // meridional arc length            (scaled by a) 
 
	// check bounds on radius of curvature
	if (r > 0.99999999*R)
		r = 0.99999999*R;
	
	if (r <= 0)
		r = 1e-12;
	
	/*-------- RESULTS FROM LUBRICATION THEORY --------*/
	
	double Dr   = R - r    ; 
	double iDr  = 1.0 /Dr  ; 
	double iDr2 = iDr*iDr  ; 
  double Dr2  = Dr  *Dr  ; 
  double Dr3  = Dr  *Dr2 ;
  double Dr4  = Dr2 *Dr2 ;
  
	double r2  = r *r ;
	double R2  = R *R ;
	double R3  = R *R2;
	double R4  = R2*R2;
	double R5  = R2*R3;
	double R6  = R3*R3;

	double Q2  = 2.0*Q*R0/R;
	double UR  = U*R;
  double g, e;
  
	// calculate local axial pressure gradient, g = (dp/dx)/Ca
	if (r > 0.0001){
  	double log = gsl_sf_log(r/R);

		g  = (8.0/(R2 - r2));
		g *= (Q2*R - U*(R2 + (R2 - r2)/(2.0*log)));
		g /= (R2 + r2 + (R2 - r2)/log);
	}
	else {
		g  = 8.0*(Q2 - UR)/R3; // limiting value as r --> 0
	}

	// calculate local membrane shear rate, e = du/dr

	// EXACT FORM:
  // double log = gsl_sf_log(r/R);
	// double rlog = r*log;
	// double irlog = 1.0/rlog;
	// e = 0.25*g*((R2 - r2)*irlog + 2.0*r) + U*irlog;
  	
	// TAYLOR SERIES ABOUT r = R:
	// (to avoid the singularity at r = 0)
	e  = iDr2 *(            - 3.0     *Q2)   ;
  e += iDr  *(2.0     *UR - 2.0     *Q2)/R ;
  e +=       (1.0     *UR - 1.45    *Q2)/R2;
  e +=  Dr  *(0.7     *UR - 1.15    *Q2)/R3;
  e +=  Dr2 *(0.55    *UR - 0.981786*Q2)/R4;
  e +=  Dr3 *(0.469286*UR - 0.882857*Q2)/R5;
  e +=  Dr4 *(0.423214*UR - 0.820839*Q2)/R6;
	// higher-order terms emphasize the singular nature of the
	// solution near the poles, so better to truncate the series
	
	/*-------------------------------------------------*/
	
	if (r > 0.0001){ // check if far from end caps
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double cosr = cos/r;
  	double sinr = sin/r;

  	// calculate function
  	f[0]  = -sin                                                                ;
		f[1]  =  cos                                                                ;
  	f[2]  = -cs                                                                 ;
  	f[3]  =  sinr*cosr + cs*sinr - qs                                           ;
  	f[4]  =  p - sig*(cs - cosr) + 0.5*(cs + cosr)*(cs*cs - cosr*cosr) + qs*sinr;
  	f[5]  =  T0*g*cos                                                           ;
  	f[6]  = -T0*e                                                               ;
  	f[7]  =  2.0*M_PI*r                                                         ;
  	f[8]  =  M_PI*r2*cos                                                        ;
		f[9]  =  (-r*ur + UR*tana)*cos/R0                                           ;
  	f[10] = -tana*cos                                                           ;
		f[11] =  M_PI*r2*x*cos/((4.0/3.0)*M_PI*v)                                   ;
		f[12] =  0.0                                                                ;
		f[13] =  0.0                                                                ;
	}
	else { // approximate cs, cphi as constant (spherical cap)
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);

  	f[0]  = -sin                                                                ;
		f[1]  =  cos                                                                ;
  	f[2]  = -cs                                                                 ;
  	f[3]  =  0.0                                                                ;
  	f[4]  =  0.5*(p - 2.0*cs*sig);                                              ;
  	f[5]  =  T0*g*cos                                                           ;
  	f[6]  = -T0*e                                                               ;
  	f[7]  =  2.0*M_PI*r                                                         ;
  	f[8]  =  M_PI*r2*cos                                                        ;
		f[9]  =  (-r*ur + UR*tana)*cos/R0                                           ;
  	f[10] = -tana*cos                                                           ;
		f[11] =  M_PI*r2*x*cos/((4.0/3.0)*M_PI*v)                                   ;
		f[12] =  0.0                                                                ;
		f[13] =  0.0                                                                ;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S;
	}
}



#endif
