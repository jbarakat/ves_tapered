/* ORDINARY DIFFERENTIAL EQUATIONS
 *  Integration schemes for a system of ODEs.
 *
 * REFERENCES
 *  Moin, Cambridge University Press (2010), pp. 64-70
 *  
 * PARAMETERS
 *  ny		[input]		number of elements in y
 *  nstep	[input]		number of integration steps
 *  par		[input]		parameters
 *  u0		[input]		source at lower boundary
 *  u1		[input]		source at upper boundary
 *  t0		[input]		lower boundary of domain
 *  t1		[input]		upper boundary of domain
 *  y0		[input]		solution vector at t = t0 (initial condition)
 *  y1		[output]	solution vector at t = t1
 */

#ifndef ODE_H
#define ODE_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "func.h"

/* PROTOTYPES */
void rk4 (int, int, double*, double, double, double, double, double *, double *);

/* IMPLEMENTATIONS */

/* Given t0 and y0 (with parameters p), integrate to y1 at t1
 * using a fourth-order Runge-Kutta scheme. */
void rk4(int ny, int nstep, double *par, double u0, double u1,
         double t0, double t1, double *y0, double *y1){
  int istep, iy, j;
  double t, u;
  double dt = (t1 - t0)/nstep;
  double du = (u1 - u0)/nstep;
  double *y00, *y01, *y02;
  double *f;

  y00 = (double*) calloc(ny , sizeof(double));
  y01 = (double*) calloc(ny , sizeof(double));
  y02 = (double*) calloc(ny , sizeof(double));
  f   = (double*) calloc(ny , sizeof(double));
  
  // initialize
  t = t0;
	u = u0;
  for (iy = 0; iy < ny; iy++)
    y02[iy] = y0[iy];

  for (istep = 0; istep < nstep; istep++){
    // setup
    t += istep*dt;
		u += istep*du;
    for (iy = 0; iy < ny; iy++)
      y00[iy] = y02[iy];

    // first step
    for (iy = 0; iy < ny; iy++)
      y01[iy] = y00[iy];
    func(ny, par, u, t, y01, f);
    for (iy = 0; iy < ny; iy++)
      y02[iy] += dt*(1.0/6.0)*f[iy];

    // second and third steps
    for (j = 0; j < 2; j++){
      for (iy = 0; iy < ny; iy++)
        y01[iy] = y00[iy] + dt*0.5*f[iy];
      func(ny, par, u, t, y01, f);
      for (iy = 0; iy < ny; iy++)
        y02[iy] += dt*(1.0/3.0)*f[iy];
    }

    // fourth step
    for (iy = 0; iy < ny; iy++)
      y01[iy] = y00[iy] + dt*f[iy];
    func(ny, par, u, t, y01, f);
    for (iy = 0; iy < ny; iy++)
      y02[iy] += dt*(1.0/6.0)*f[iy];
  }

  // store results
  for (iy = 0; iy < ny; iy++)
    y1[iy] = y02[iy];
  
  free(y00);
  free(y01);
  free(y02);
  free(f  );
}

#endif
