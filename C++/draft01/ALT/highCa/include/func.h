/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  par					[input]		parameters
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

/* IMPLEMENTATIONS */
void func(int ny, double *par,
          double t, double *y, double *f){
  int i;
  if (ny != 12)
    cout << "Error: ny should equal 12." << endl;

	// define parameters
	double v    = par[0];
	double R0   = par[1];
	double alph = par[2];
	double ur   = 0.0;

	double area = 4.0*M_PI;
	double vlme = 4.0*M_PI*v/3.0;
	double tana = gsl_sf_sin(alph)/gsl_sf_cos(alph);
  
	// define variables
  double r   = y[0 ];
  double psi = y[1 ];
  double p   = y[2 ];
  double sig = y[3 ];
  double A   = y[4 ];
  double V   = y[5 ];
  double Q2  = y[6 ];
  double R   = y[7 ];
	double U   = y[8 ];
  double S   = y[9 ];
	double x   = y[10];
	double xcm = y[11];
  
	double R2 = R*R;
	double UR = U*R;
	
	// check bounds on radius
	if (r > 0.999999999*R)
		r = 0.999999999*R;
	
	if (r <= 0)
		r = 1e-12;
	
	if (r > 0.002){ // check if far from end caps
  	double r2  = r*r;
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double log = gsl_sf_log(r/R);
  	double cosr = cos/r;
  	double sinr = sin/r;
  	double g, e;

  	// calculate local flow coefficient
  	g  = (8.0/(R2 - r2));
		g *= (Q2*R - U*(R2 + (R2 - r2)/(2.0*log)));
  	g /= (R2 + r2 + (R2 - r2)/log);

  	// calculate local shear rate
  	// (expand e in Taylor series about r = 1)
  	double Dr  = R - r; 
  	double Dr2 = Dr*Dr; 
  	double Dr3 = Dr2*Dr;
  	double Dr4 = Dr2*Dr2;

  	e  = -(3.0/Dr2)*Q2;
  	e += (2.0/Dr )*(UR - Q2);
  	e += UR - (29.0/20.0)*Q2;
  	e += Dr *(7.0/10.0*UR - (32.0/20.0)*Q2);
  	e += Dr2*(11.0/20.0*UR - (2749.0/2800.0)*Q2);
  	e += Dr3*(657.0/1400.0*UR - (309.0/350.0)*Q2);
  	e += Dr4*(237.0/560.0*UR - (45967.0/56000.0)*Q2);

  	// calculate function
  	f[0 ] = -sin;
  	f[1 ] = -p/sig - cosr;
  	f[2 ] = g*cos;
  	f[3 ] = -e;
  	f[4 ] = 2.0*M_PI*r;
  	f[5 ] = M_PI*r2*cos;
		f[6 ] = (-r*ur + UR*tana)*cos/R0;
  	f[7 ] = -tana*cos;
  	f[8 ] = 0.0;
  	f[9 ] = 0.0;
		f[10] = cos;
  	f[11] = M_PI*r*r*x*cos/vlme;
	}
	else { // approximate p, sig as constant 
	       // and set dpsi/ds = -cos(psi)/r
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
		
		f[0 ] = -sin;
		f[1 ] = -p/(2.0*sig);
		f[2 ] = 0.0;
		f[3 ] = 0.0;
		f[4 ] = 2.0*M_PI*r;
		f[5 ] = M_PI*r*r*cos;
		f[6 ] = (-r*ur + UR*tana)*cos/R0;
  	f[7 ] = -tana*cos;
  	f[8 ] = 0.0;
  	f[9 ] = 0.0;
		f[10] = cos;
  	f[11] = M_PI*r*r*x*cos/vlme;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S;
	}
}



#endif
