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
#include <gsl/gsl_sf_trig.h>

using namespace std;

/* PROTOTYPES */
/* Note: input and output files have different functions because
 * of the naming convention:
 * - input files labeled by CaXXX.dat suffix
 * - output files all use CaInf.dat suffix */
void fileCheckInput(int, int, double, bool &); 
void fileCheckOutput(int, int, bool &); 
void readInput(int, int, int, int, double, double *, double *); 
void readOutput(int, int, int, int, double *, double *); 

/* IMPLEMENTATIONS */
/* Check if input file exists. */
void fileCheckInput(int v, int conf, double Ca, bool &info){
  string filename, dir;
	
	// convert double to string
	ostringstream capnum;
	if (Ca >= 0.000000001 && Ca < 0.00000001)
		capnum << fixed << setprecision(10) << Ca;
	if (Ca >= 0.00000001 && Ca < 0.0000001)
		capnum << fixed << setprecision(9) << Ca;
	if (Ca >= 0.0000001 && Ca < 0.000001)
		capnum << fixed << setprecision(8) << Ca;
	if (Ca >= 0.000001 && Ca < 0.00001)
		capnum << fixed << setprecision(7) << Ca;
	if (Ca >= 0.00001 && Ca < 0.0001)
		capnum << fixed << setprecision(6) << Ca;
	if (Ca >= 0.0001 && Ca < 0.001)
		capnum << fixed << setprecision(5) << Ca;
	if (Ca >= 0.001 && Ca < 0.01)
		capnum << fixed << setprecision(4) << Ca;
	if (Ca >= 0.01 && Ca < 0.1)
		capnum << fixed << setprecision(3) << Ca;
	if (Ca >= 0.1 && Ca < 1.0)
		capnum << fixed << setprecision(2) << Ca;
	if (Ca >= 1.0)
		capnum << fixed << setprecision(0) << Ca;

  // get filename
  ostringstream redvol, confin;
  redvol << fixed << v;
  confin << fixed << conf;
  dir = "../input";
  filename = "/v" + redvol.str() + "/conf" + confin.str()
             + "/sln_v" + redvol.str() + "_conf" + confin.str()
             + "_Ca" + capnum.str();
  filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
  filename = dir + filename + ".dat";

  // check if file exists
  ifstream file(filename.c_str());
  if (file)
    info = true; // file exists
  else
    info = false; // file does not exist
}

/* Check if output file exists. */
void fileCheckOutput(int v, int conf, bool &info){
  string filename, dir;

  // get filename
  ostringstream redvol, confin;
  redvol << fixed << v;
  confin << fixed << conf;
  dir = "../output";
  filename = "/v" + redvol.str() + "/conf" + confin.str()
             + "/sln_v" + redvol.str() + "_conf" + confin.str()
             + "_CaInf";
  filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
  filename = dir + filename + ".dat";

  // check if file exists
  ifstream file(filename.c_str());
  if (file)
    info = true; // file exists
  else
    info = false; // file does not exist
}


/* Read input .dat file containing the solution s. */
void readInput(int n, int m, int v, int conf, double Ca, double *t, double *s){
  // error flags
  if (n != 11){
    cout << "Error: support only for n = 11." << endl;
		return;
  }

//  // error flags
//  if (n != 11){
//    cout << "Error: support only for n = 11." << endl;
//		return;
//  }

 // if (v != 65 && v != 70 && v != 75 && v!=80 && v != 85 && v != 90 && v != 95){
 //   cout << "Error: support only for v = 85, 90, 95." << endl;
 // 	return;
 // }

 // if (conf < 90 || conf > 99){
 //   cout << "Error: support only for 90 <= conf <= 99." << endl;
 // 	return;
 // }

  typedef vector<double> vdouble;

  // declare variables
  int    i, j, k, i1, j1, k1; 
  string filename, line, dir;
  ifstream file;
  vdouble T;
  vdouble R   , X   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , UU  , RRt  , XXcm  , SS   ;
  vdouble r(m), x(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), U(m), Rt(m), xcm(m), S (m);
	
	// convert double to string
	ostringstream capnum;
	if (Ca >= 0.000000001 && Ca < 0.00000001)
		capnum << fixed << setprecision(10) << Ca;
	if (Ca >= 0.00000001 && Ca < 0.0000001)
		capnum << fixed << setprecision(9) << Ca;
	if (Ca >= 0.0000001 && Ca < 0.000001)
		capnum << fixed << setprecision(8) << Ca;
	if (Ca >= 0.000001 && Ca < 0.00001)
		capnum << fixed << setprecision(7) << Ca;
	if (Ca >= 0.00001 && Ca < 0.0001)
		capnum << fixed << setprecision(6) << Ca;
	if (Ca >= 0.0001 && Ca < 0.001)
		capnum << fixed << setprecision(5) << Ca;
	if (Ca >= 0.001 && Ca < 0.01)
		capnum << fixed << setprecision(4) << Ca;
	if (Ca >= 0.01 && Ca < 0.1)
		capnum << fixed << setprecision(3) << Ca;
	if (Ca >= 0.1 && Ca < 1.0)
		capnum << fixed << setprecision(2) << Ca;
	if (Ca >= 1.0)
		capnum << fixed << setprecision(0) << Ca;

  // get filename
  ostringstream redvol, confin;
  redvol << fixed << v;
  confin << fixed << conf;
  dir = "../input";
  filename = "/v" + redvol.str() + "/conf" + confin.str()
             + "/sln_v" + redvol.str() + "_conf" + confin.str()
             + "_Ca" + capnum.str();
  filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
  filename = dir + filename + ".dat";
  cout << "Reading " << filename << "." << endl;

	// open file
  file.open(filename.c_str());
  while (file){
    getline(file, line);
    istringstream iss(line);
    vdouble row;

    for (i = 0; i < 11; i++){
      double data;
      iss >> data;
      row.push_back(data);
    }

    if (file.eof())
      break;

    T  .push_back(row[0 ]);
    R  .push_back(row[1 ]);
    PSI.push_back(row[2 ]);
    CS .push_back(row[3 ]);
    QS .push_back(row[4 ]);
    P  .push_back(row[5 ]);
    SIG.push_back(row[6 ]);
    AA .push_back(row[7 ]);
    VV .push_back(row[8 ]);
    QQ .push_back(row[9 ]);
    SS .push_back(row[10]);
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
    r  [i] = R  [k] + (R  [k1] - R  [k])*(Dt/DT);
    psi[i] = PSI[k] + (PSI[k1] - PSI[k])*(Dt/DT);
    cs [i] = CS [k] + (CS [k1] - CS [k])*(Dt/DT);
    qs [i] = QS [k] + (QS [k1] - QS [k])*(Dt/DT);
    p  [i] = P  [k] + (P  [k1] - P  [k])*(Dt/DT);
    sig[i] = SIG[k] + (SIG[k1] - SIG[k])*(Dt/DT);
    A  [i] = AA [k] + (AA [k1] - AA [k])*(Dt/DT);
    V  [i] = VV [k] + (VV [k1] - VV [k])*(Dt/DT);
    Q  [i] = QQ [k] + (QQ [k1] - QQ [k])*(Dt/DT);
		Rt [i] = 1.0                                ;
		U  [i] = 1.0                                ;
    S  [i] = SS [k] + (SS [k1] - SS [k])*(Dt/DT);
  }

	// calculate x and xcm using by forward differences
	double ds, dx, dxcm;
	x  [0] = 0.0;
	xcm[0] = 0.0;
	for (i = 0; i < m-1; i++){
		i1 = i + 1;
		ds   = S[0]*(t[i1] - t[i]);
		dx   = gsl_sf_cos(psi[i])*ds;
		dxcm = M_PI*r[i]*r[i]*x[i]*dx/V[m-1];

		x  [i1] = x  [i] + dx  ;
		xcm[i1] = xcm[i] + dxcm;
	}

//	// shift origin to center of mass, then recalculate x and xcm
//	for (j = 0; j < 5; j++){ // repeat several times to improve accuracy
//		for (i = 0; i < m; i++){
//			x[i] = x[i] - xcm[m-1];
//		}
//		
//		for (i = 0; i < m-1; i++){
//			i1 = i + 1;
//			ds   = S[0]*(t[i1] - t[i]);
//			dx   = gsl_sf_cos(psi[i])*ds;
//			dxcm = M_PI*r[i]*r[i]*x[i]*dx/V[m-1];
//	
//			x  [i1] = x  [i] + dx  ;
//			xcm[i1] = xcm[i] + dxcm;
//		}
//	}
  
	// assemble the solution vector
  for (i = 0; i < m; i++){
    s[i*n + 0 ] = r  [i];
//    s[i*n + 1 ] = x  [i];
    s[i*n + 1 ] = psi[i];
    s[i*n + 2 ] = p  [i];
    s[i*n + 3 ] = sig[i];
    s[i*n + 4 ] = A  [i];
    s[i*n + 5 ] = V  [i];
    s[i*n + 6 ] = Q  [i];
    s[i*n + 7 ] = Rt [i];
//    s[i*n + 9 ] = xcm[i];
    s[i*n + 8 ] = U  [i];
    s[i*n + 9 ] = S  [i];
    s[i*n + 10] = x  [i];
  }

//  // assemble the solution vector
//  for (i = 0; i < m; i++){
//    s[i*n + 0 ] = r  [i];
////    s[i*n + 1 ] = x  [i];
//    s[i*n + 2 ] = psi[i];
//    s[i*n + 3 ] = p  [i];
//    s[i*n + 4 ] = sig[i];
//    s[i*n + 5 ] = A  [i];
//    s[i*n + 6 ] = V  [i];
//    s[i*n + 7 ] = Q  [i];
//    s[i*n + 8 ] = Rt [i];
////    s[i*n + 9 ] = xcm[i];
//    s[i*n + 9] = U  [i];
//    s[i*n + 10] = S  [i];
//  }
}


/* Read output .dat file containing the solution s. */
void readOutput(int n, int m, int v, int conf, double *t, double *s){
  // error flags
  if (n != 8){
    cout << "Error: support only for n = 8." << endl;
		return;
  }

 // if (v != 65 && v != 70 && v != 75 && v!=80 && v != 85 && v != 90 && v != 95){
 //   cout << "Error: support only for v = 85, 90, 95." << endl;
 // 	return;
 // }

 // if (conf < 90 || conf > 99){
 //   cout << "Error: support only for 90 <= conf <= 99." << endl;
 // 	return;
 // }

  typedef vector<double> vdouble;

  // declare variables
  int    i, j, k, i1, j1, k1; 
  string filename, line, dir;
  ifstream file;
  vdouble T;
  vdouble R   , PSI   , CS   , QS   , P   , SIG   , AA   , VV   , QQ   , SS   ;   
  vdouble r(m), psi(m), cs(m), qs(m), p(m), sig(m), A (m), V (m), Q (m), S (m);
	
  // get filename
  ostringstream redvol, confin;
  redvol << fixed << v;
  confin << fixed << conf;
  dir = "../output";
  filename = "/v" + redvol.str() + "/conf" + confin.str()
             + "/sln_v" + redvol.str() + "_conf" + confin.str()
             + "_CaInf";
  filename.erase(remove(filename.begin(), filename.end(), '.'), filename.end());
  filename = dir + filename + ".dat";
  cout << "Reading " << filename << "." << endl;

	// open file
  file.open(filename.c_str());
  while (file){
    getline(file, line);
    istringstream iss(line);
    vdouble row;

    for (i = 0; i < 11; i++){
      double data;
      iss >> data;
      row.push_back(data);
    }

    if (file.eof())
      break;

    T  .push_back(row[0 ]);
    R  .push_back(row[1 ]);
    PSI.push_back(row[2 ]);
    CS .push_back(row[3 ]);
    QS .push_back(row[4 ]);
    P  .push_back(row[5 ]);
    SIG.push_back(row[6 ]);
    AA .push_back(row[7 ]);
    VV .push_back(row[8 ]);
    QQ .push_back(row[9 ]);
    SS .push_back(row[10]);
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
    r  [i] = R  [k] + (R  [k1] - R  [k])*(Dt/DT);
    psi[i] = PSI[k] + (PSI[k1] - PSI[k])*(Dt/DT);
    cs [i] = CS [k] + (CS [k1] - CS [k])*(Dt/DT);
    qs [i] = QS [k] + (QS [k1] - QS [k])*(Dt/DT);
    p  [i] = P  [k] + (P  [k1] - P  [k])*(Dt/DT);
    sig[i] = SIG[k] + (SIG[k1] - SIG[k])*(Dt/DT);
    A  [i] = AA [k] + (AA [k1] - AA [k])*(Dt/DT);
    V  [i] = VV [k] + (VV [k1] - VV [k])*(Dt/DT);
    Q  [i] = QQ [k] + (QQ [k1] - QQ [k])*(Dt/DT);
    S  [i] = SS [k] + (SS [k1] - SS [k])*(Dt/DT);
  }

  // assemble the solution vector
  for (i = 0; i < m; i++){
    s[i*n + 0] = r  [i];
    s[i*n + 1] = psi[i];
    s[i*n + 2] = p  [i];
    s[i*n + 3] = sig[i];
    s[i*n + 4] = A  [i];
    s[i*n + 5] = V  [i];
    s[i*n + 6] = Q  [i];
    s[i*n + 7] = S  [i];
  }
}



#endif
