/* WRITE FILE
 *  Write .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  m			[input]		number of shooting points
 *  id		[input]		identifier for writing files
 *  t			[output]	abscissas
 *  s			[input]		solution vector
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
void writeSoln(int n, int m, int *id,
               double *t, double *s, string path){
	// error flag
	if (n != 12){
		cout << "Error: support only for n = 12." << endl;
		return;
	}

	// declare variables
	int i;
	int width = 12;
	ostringstream ssRedVol, ssTapAng;
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
	int v    = id[0];
	int alph = id[1];
	ssRedVol << v   ;
	ssTapAng << alph;
	
	// set directory
	dir = path + "/v"    + ssRedVol.str() 
	           + "/a00"  + ssTapAng.str();
	
	// set filename
	fn = "/sln_v" + ssRedVol.str()
	   + "_a00"  + ssTapAng.str() 
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
	       << setw(width)   << "p" 
				 << setw(width)   << "sig"
				 << setw(width)   << "A"
				 << setw(width)   << "V"
				 << setw(width)   << "Q2"
				 << setw(width)   << "R"
				 << setw(width)   << "U"
				 << setw(width)   << "S"
				 << setw(width)   << "x"
				 << setw(width)   << "xcm";
	file << header.str() << endl;
	for (i = 0; i < m; i++){
		ostringstream tt, r, psi, cs, qs, p, 
		                  sig, A, V , Q , S,
											U, R, x, xcm;
		tt  << fixed << setw(width)   << t[i];
		r   << fixed << setw(width)   << s[i*n + 0 ];
		psi << fixed << setw(width)   << s[i*n + 1 ];
		cs  << fixed << setw(width)   << CS[i]      ;
		qs  << fixed << setw(width)   << QS[i]      ;
		p   << fixed << setw(width)   << s[i*n + 2 ];
		sig << fixed << setw(width)   << s[i*n + 3 ];
		A   << fixed << setw(width)   << s[i*n + 4 ];
		V   << fixed << setw(width)   << s[i*n + 5 ];
		Q   << fixed << setw(width)   << s[i*n + 6 ];
		R   << fixed << setw(width)   << s[i*n + 7 ];
		U   << fixed << setw(width)   << s[i*n + 8 ];
		S   << fixed << setw(width)   << s[i*n + 9 ];
		x   << fixed << setw(width)   << s[i*n + 10];
		xcm << fixed << setw(width)   << s[i*n + 11];
		
		line = tt.str() + r.str() + psi.str() + cs.str() + qs.str()
		     + p.str() + sig.str() + A.str() + V .str() + Q .str() 
				 + R.str() + U.str() + S.str() + x.str() + xcm.str();

		file << line << endl;
	}
	file.close();

	cout << "Solution written to " << fn << "." << endl;
}


#endif
