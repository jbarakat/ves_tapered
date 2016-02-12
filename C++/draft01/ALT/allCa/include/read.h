/* READ FILE
 *  Read .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  Seifert et al, Phys Rev A (1991)
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

#ifndef READ_H
#define READ_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */
void fileCheck(int, int, double, bool &);
void readInput(int, int, int *, double *, double *);
void readOutput(int, int, int, int, double, double *, double *);

/* IMPLEMENTATIONS */

/* Check if output file exists. */
//void fileCheck(int v, int conf, double Ca, bool &info){
//	string filename, dir;
//	
//	// convert double to string
//	ostringstream capnum;
//	if (Ca >= 0.000000001 && Ca < 0.00000001)
//		capnum << fixed << setprecision(10) << Ca;
//	if (Ca >= 0.00000001 && Ca < 0.0000001)
//		capnum << fixed << setprecision(9) << Ca;
//	if (Ca >= 0.0000001 && Ca < 0.000001)
//		capnum << fixed << setprecision(8) << Ca;
//	if (Ca >= 0.000001 && Ca < 0.00001)
//		capnum << fixed << setprecision(7) << Ca;
//	if (Ca >= 0.00001 && Ca < 0.0001)
//		capnum << fixed << setprecision(6) << Ca;
//	if (Ca >= 0.0001 && Ca < 0.001)
//		capnum << fixed << setprecision(5) << Ca;
//	if (Ca >= 0.001 && Ca < 0.01)
//		capnum << fixed << setprecision(4) << Ca;
//	if (Ca >= 0.01 && Ca < 0.1)
//		capnum << fixed << setprecision(3) << Ca;
//	if (Ca >= 0.1 && Ca < 1.0)
//		capnum << fixed << setprecision(2) << Ca;
//	if (Ca >= 1.0)
//		capnum << fixed << setprecision(0) << Ca;
// 
//	// get filename
//	ostringstream redvol, confin;
//	redvol << fixed << v;
//	confin << fixed << conf;
//	dir = "../output";
//	filename = "/v" + redvol.str() + "/conf" + confin.str()
//	           + "/sln_v" + redvol.str() + "_conf" + confin.str()
//	           + "_Ca" + capnum.str();
//	filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
//	filename = dir + filename + ".dat";
//	
//	// check if file exists
//	ifstream file(filename.c_str());
//	if (file)
//		info = true; // file exists
//	else
//		info = false; // file does not exist
//}

/* Read input .dat file (e.g. generated from straight-channel problem). */
void readInput(int n, int m, int *id, double *t, double *s){
	// error flags
	if (n != 14){
		cout << "Error: support only for n = 14." << endl;
		return;
	}

	int v    = id[0];
	int conf = id[1];

	typedef vector<double> vdouble;

	// declare variables
	int    i, j, k, i1, j1, k1;
	string filename, line;
	ifstream file;
	vdouble T;
	vdouble R   , X   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , RTUBE   , XCOM   , UU  , SS   ;
	vdouble r(m), x(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), rtube(m), xcom(m), U(m), S (m);
 
	// get filename
	ostringstream redvol, confin;
	redvol << fixed << v;
	confin << fixed << conf;
	filename = "../input/sln_v" + redvol.str() + "_conf" + confin.str() + "_kb0.dat";
	cout << "Reading " << filename << "." << endl;

	// open file
	file.open(filename.c_str());
	while (file){
		int skiplines = 1;
		for (i = 0; i < skiplines; i++)
			getline(file, line);
		istringstream iss(line);
		vdouble row;
		
		for (i = 0; i < n+1; i++){
			double data;
			iss >> data;
			row.push_back(data);
		}
		
		if (file.eof())
			break;

		T    .push_back(row[0 ]);
		R    .push_back(row[1 ]);
		X    .push_back(row[2 ]);
		PSI  .push_back(row[3 ]);
		CS   .push_back(row[4 ]);
		QS   .push_back(row[5 ]);
		P    .push_back(row[6 ]);
		SIG  .push_back(row[7 ]);
		AA   .push_back(row[8 ]);
		VV   .push_back(row[9 ]);
		QQ   .push_back(row[10]);
		RTUBE.push_back(row[11]);
		XCOM .push_back(row[12]);
		UU   .push_back(row[13]);
		SS   .push_back(row[14]);
	}

	file.close();

	// interpolate
	int    nt  = T.size() - 1;
	double t0 = 0.0;
	double t1 = 1.0;
	double dt = (t1 - t0)/(m-1);
	double ti;
	double DT, Dt;

	// assemble the abscissas
	for (i = 0; i < m; i++){
		t[i] = i*dt;
	}

	for (i = 0; i < m; i++){
		ti = t[i];

		// find ti in T
		for (j = 0; j < nt; j++){
			j1 = j + 1;
			if (ti >= T[j] && ti < T[j1] ){
				k = j;
				break;
			}
			else if (ti == T[n]){
				k = nt-1;
				break;
			}
		}

		k1 = k + 1;
		DT = T[k1] - T[k];
		Dt = ti    - T[k];

		// linear interpolation
		r    [i] = R    [k] + (R    [k1] - R    [k])*(Dt/DT);
		x    [i] = X    [k] + (X    [k1] - X    [k])*(Dt/DT);
		psi  [i] = PSI  [k] + (PSI  [k1] - PSI  [k])*(Dt/DT);
		cs   [i] = CS   [k] + (CS   [k1] - CS   [k])*(Dt/DT);
		qs   [i] = QS   [k] + (QS   [k1] - QS   [k])*(Dt/DT);
		p    [i] = P    [k] + (P    [k1] - P    [k])*(Dt/DT);
		sig  [i] = SIG  [k] + (SIG  [k1] - SIG  [k])*(Dt/DT);
		A    [i] = AA   [k] + (AA   [k1] - AA   [k])*(Dt/DT);
		V    [i] = VV   [k] + (VV   [k1] - VV   [k])*(Dt/DT);
		Q    [i] = QQ   [k] + (QQ   [k1] - QQ   [k])*(Dt/DT);
		rtube[i] = RTUBE[k] + (RTUBE[k1] - RTUBE[k])*(Dt/DT);
		xcom [i] = XCOM [k] + (XCOM [k1] - XCOM [k])*(Dt/DT);
		U    [i] = UU   [k] + (UU   [k1] - UU   [k])*(Dt/DT);
		S    [i] = SS   [k] + (SS   [k1] - SS   [k])*(Dt/DT);
	}

	// assemble the solution vector
	for (i = 0; i < m; i++){
    s[i*n + 0 ] = r    [i];
    s[i*n + 1 ] = x    [i];
    s[i*n + 2 ] = psi  [i];
    s[i*n + 3 ] = cs   [i];
    s[i*n + 4 ] = qs   [i];
    s[i*n + 5 ] = p    [i];
    s[i*n + 6 ] = sig  [i];
    s[i*n + 7 ] = A    [i];
    s[i*n + 8 ] = V    [i];
    s[i*n + 9 ] = Q    [i];
    s[i*n + 10] = rtube[i];		
    s[i*n + 11] = xcom [i];		
    s[i*n + 12] = U    [i];		
    s[i*n + 13] = S    [i];		
	}
}

///* Read output .dat file containing the solution s. */
//void readOutput(int n, int m, int v, int conf, double Ca, double *t, double *s){
//	// error flags
//	if (n != 10){
//		cout << "Error: support only for n = 10." << endl;
//		return;
//	}
//
//	typedef vector<double> vdouble;
//
//	// declare variables
//	int    i, j, k, i1, j1, k1;
//	string filename, line, dir;
//	ifstream file;
//	vdouble T;
//	vdouble R   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , SS   ;
//	vdouble r(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), S (m);
//	
//	// convert double to string
//	ostringstream capnum;
//	if (Ca >= 0.000000001 && Ca < 0.00000001)
//		capnum << fixed << setprecision(10) << Ca;
//	if (Ca >= 0.00000001 && Ca < 0.0000001)
//		capnum << fixed << setprecision(9) << Ca;
//	if (Ca >= 0.0000001 && Ca < 0.000001)
//		capnum << fixed << setprecision(8) << Ca;
//	if (Ca >= 0.000001 && Ca < 0.00001)
//		capnum << fixed << setprecision(7) << Ca;
//	if (Ca >= 0.00001 && Ca < 0.0001)
//		capnum << fixed << setprecision(6) << Ca;
//	if (Ca >= 0.0001 && Ca < 0.001)
//		capnum << fixed << setprecision(5) << Ca;
//	if (Ca >= 0.001 && Ca < 0.01)
//		capnum << fixed << setprecision(4) << Ca;
//	if (Ca >= 0.01 && Ca < 0.1)
//		capnum << fixed << setprecision(3) << Ca;
//	if (Ca >= 0.1 && Ca < 1.0)
//		capnum << fixed << setprecision(2) << Ca;
//	if (Ca >= 1.0)
//		capnum << fixed << setprecision(0) << Ca;
// 
//	// get filename
//	ostringstream redvol, confin;
//	redvol << fixed << v;
//	confin << fixed << conf;
//	dir = "../output";
//	filename = "/v" + redvol.str() + "/conf" + confin.str()
//	           + "/sln_v" + redvol.str() + "_conf" + confin.str()
//	           + "_Ca" + capnum.str();
//	filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
//	filename = dir + filename + ".dat";
//	cout << "Reading " << filename << "." << endl;
//
//	// open file
//	file.open(filename.c_str());
//	while (file){
//		getline(file, line);
//		istringstream iss(line);
//		vdouble row;
//		
//		for (i = 0; i < 11; i++){
//			double data;
//			iss >> data;
//			row.push_back(data);
//		}
//		
//		if (file.eof())
//			break;
//
//		T  .push_back(row[0 ]);
//		R  .push_back(row[1 ]);
//		PSI.push_back(row[2 ]);
//		CS .push_back(row[3 ]);
//		QS .push_back(row[4 ]);
//		P  .push_back(row[5 ]);
//		SIG.push_back(row[6 ]);
//		AA .push_back(row[7 ]);
//		VV .push_back(row[8 ]);
//		QQ .push_back(row[9 ]);
//		SS .push_back(row[10]);
//	}
//
//	file.close();
//
//	// interpolate
//	int    nt  = T.size() - 1;
//	double t0 = 0.0;
//	double t1 = 1.0;
//	double dt = (t1 - t0)/(m-1);
//	double ti;
//	double DT, Dt;
//
//	// assemble the abscissas
//	for (i = 0; i < m; i++){
//		t[i] = i*dt;
//	}
//
//	for (i = 0; i < m; i++){
//		ti = t[i];
//
//		// find ti in T
//		for (j = 0; j < nt; j++){
//			j1 = j + 1;
//			if (ti >= T[j] && ti < T[j1] ){
//				k = j;
//				break;
//			}
//			else if (ti == T[n]){
//				k = nt-1;
//				break;
//			}
//		}
//
//		k1 = k + 1;
//		DT = T[k1] - T[k];
//		Dt = ti    - T[k];
//
//		// linear interpolation
//		r  [i] = R  [k] + (R  [k1] - R  [k])*(Dt/DT);
//		psi[i] = PSI[k] + (PSI[k1] - PSI[k])*(Dt/DT);
//		cs [i] = CS [k] + (CS [k1] - CS [k])*(Dt/DT);
//		qs [i] = QS [k] + (QS [k1] - QS [k])*(Dt/DT);
//		p  [i] = P  [k] + (P  [k1] - P  [k])*(Dt/DT);
//		sig[i] = SIG[k] + (SIG[k1] - SIG[k])*(Dt/DT);
//		A  [i] = AA [k] + (AA [k1] - AA [k])*(Dt/DT);
//		V  [i] = VV [k] + (VV [k1] - VV [k])*(Dt/DT);
//		Q  [i] = QQ [k] + (QQ [k1] - QQ [k])*(Dt/DT);
//		S  [i] = SS [k] + (SS [k1] - SS [k])*(Dt/DT);
//	}
//
//	// assemble the solution vector
//	for (i = 0; i < m; i++){
//    s[i*n + 0] = r  [i];
//    s[i*n + 1] = psi[i];
//    s[i*n + 2] = cs [i];
//    s[i*n + 3] = qs [i];
//    s[i*n + 4] = p  [i];
//    s[i*n + 5] = sig[i];
//    s[i*n + 6] = A  [i];
//    s[i*n + 7] = V  [i];
//    s[i*n + 8] = Q  [i];
//    s[i*n + 9] = S  [i];		
//	}
//}

#endif
