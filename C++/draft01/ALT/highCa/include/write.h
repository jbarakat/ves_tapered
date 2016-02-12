/* WRITE FILE
 *  Write .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  m			[input]		number of shooting points
 *  v			[input]		reduced volume (x 100)
 *  conf	[input]		confinement parameter scaled by critical value (x 100)
 *  t			[output]	abscissas (meridional arc length scaled to the [0,1] axis)
 *  r			[output]	cylindrical radius
 *  psi		[output]	tilt angle
 *  cs		[output]	meridional curvature (points outward from closed contour)
 *  qs		[output]	meridional component of transverse shear tension
 *  p			[output]	pressure difference between exterior and interior
 *  sig		[output]	mean tension
 *  A			[output]	surface area
 *  V			[output]	volume
 *	Q			[output]	leakback flow rate
 *  S			[output]	total meridional arc length (half-space)
 */

#ifndef WRITE_H
#define WRITE_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_sf_trig.h>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Write .dat file. */
void writeSoln(int n, int m, int v, int conf,
               double *t, double *s, string path){
	// error flag
	if (n != 8){
		cout << "Error: support only for n = 8." << endl;
		return;
	}

	// declare variables
	int i;
	int width = 14;
	ostringstream ssRedVol, ssConfin;
	ostringstream header;
	ofstream file;
	string dir, fn, line;

	// set transverse shear tension to zero everywhere
	vector<double> QS(m, 0.0);

	// calculate meridional curvature using Laplace relation
	vector<double> CS(m);
	double rr, ppsi, pp, ssig;
	for (i = 0; i < m; i++){
		rr   = s[i*n + 0];
		ppsi = s[i*n + 1];
		pp   = s[i*n + 2];
		ssig = s[i*n + 3];
		if (rr > 0.02)
			CS[i] = pp/ssig + gsl_sf_cos(ppsi)/rr;
		else
			CS[i] = pp/(2.0*ssig);
	}

	// convert integers to strings
	ssRedVol << v   ;
	ssConfin << conf;
	
	// set directory
	dir = path + "/v"    + ssRedVol.str() 
	           + "/conf" + ssConfin.str();
	
	// set filename
	fn = "/sln_v" + ssRedVol.str()
	   + "_conf"  + ssConfin.str() 
		 + "_CaInf" ;
	fn.erase(remove(fn.begin(), fn.end(), '.'), fn.end());
	fn = dir + fn + ".dat";

	// write to file
	file.open(fn.c_str());
	header << setw(width)   << "t" 
	       << setw(width)   << "r" 
				 << setw(width)   << "psi"
				 << setw(width)   << "cs"
				 << setw(width)   << "qs"
	       << setw(width+4) << "p" 
				 << setw(width+4) << "sig"
				 << setw(width)   << "A"
				 << setw(width)   << "V"
				 << setw(width)   << "Q2"
				 << setw(width)   << "S";
	file << header.str() << endl;
	for (i = 0; i < m; i++){
		ostringstream tt, r, psi, cs, qs, p, 
		                  sig, A, V , Q , S;
		tt  << fixed << setw(width)   << t[i];
		r   << fixed << setw(width)   << s[i*n + 0];
		psi << fixed << setw(width)   << s[i*n + 1];
		cs  << fixed << setw(width)   << CS[i]     ;
		qs  << fixed << setw(width)   << QS[i]     ;
		p   << fixed << setw(width+4) << s[i*n + 2];
		sig << fixed << setw(width+4) << s[i*n + 3];
		A   << fixed << setw(width)   << s[i*n + 4];
		V   << fixed << setw(width)   << s[i*n + 5];
		Q   << fixed << setw(width)   << s[i*n + 6];
		S   << fixed << setw(width)   << s[i*n + 7];
		
		line = tt.str() + r.str() + psi.str() + cs.str() + qs.str()
		     + p.str() + sig.str() + A.str() + V .str() + Q .str() + S.str();

		file << line << endl;
	}
	file.close();

	cout << "Solution written to " << fn << "." << endl;
}


#endif
