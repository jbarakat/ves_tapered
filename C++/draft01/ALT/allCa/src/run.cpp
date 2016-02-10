/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include "../include/shoot.h"
#include "../include/read.h"
#include "../include/write.h"

void init(int, int, int, int, double, double&, double&, double*, double*);
void getCa(bool, vector<double>&);
void getVr(bool, vector<int>&);
void getCf(bool, vector<int>&);
void getCritCf(int, double&);

int main(){
	int    i, j, k, l;
	int    n    = 10;			// size of solution vector
	int    m    = 403;		// number of shooting points
	int    maxiter = 251;	// maximum number of iterations
	double tol = 1e-8;		// tolerance
	int    nrk;						// number of Runge-Kutta steps

	int    v;							// reduced volume
	int    conf;					// confinement
	int    flag;					// error flag
	bool   info;					// boolean for file check

	// output directory
	string opath = "../output";

	// abscissa and solution vectors
	vector<double> t(m), si(n*m), sm(n*m), sf(n*m);
	
	// parameters
	double Ca, area, vlme;
	double dC0, dC1, slope;
	vector<double> vecCa;
	vector<int   > vecVr;
	vector<int   > vecCf;

	// get vector of capillary numbers
	bool flipCa = true;
	getCa(flipCa, vecCa);

	// get vector of reduced volumes
	bool flipVr = false;
	getVr(flipVr, vecVr);

	// get vector of confinement parameters
	bool flipCf = false;
	getCf(flipCf, vecCf);
	
	// loop over reduced volume
	for (l = 0; l < vecVr.size(); l++){
		v = vecVr[l];

		// loop over confinement parameter
		for (k = 0; k < vecCf.size(); k++){
			conf = vecCf[k];
			
			// initialize
			nrk = 60;
			cout << "Initializing... " << endl;
			init(n, m, v, conf, vecCa[0], area, vlme, t.data(), si.data());
			for (j = 0; j < m*n; j++){
				sm[j] = si[j];
			}
			cout << "Initialization complete." << endl;
			
			// loop over capillary numbers
			for (i = 0; i < vecCa.size(); i++){
				if (Ca > 70)
					nrk = 100;

				// choose capillary number
				Ca = vecCa[i];

				// check if file exists
				fileCheck(v, conf, vecCa[i], info);
				if (info){ // file exists
					// update solution
					readOutput(n, m, v, conf, vecCa[i], t.data(), si.data());
					for (j = 0; j < m*n; j++){
						sm[j] = si[j];
					}
			
			//		// skip to next step
			//		cout << "Output file already exists. Skipping..." << endl;
			//		continue;
				}
				else { // check neighboring solutions
			//		// or just skip to next step
			//		cout << "Output file does not exist. Skipping..." << endl;
			//		continue;

					bool info1, info2, info3, info4, info5, info6;
					if (i == 0)
						fileCheck(v, conf, vecCa[i], info1);
					else
						fileCheck(v, conf, vecCa[i-1], info1);
					if (i == vecCa.size()-1)
						fileCheck(v, conf, vecCa[i], info2);
					else
						fileCheck(v, conf, vecCa[i+1], info2);
					fileCheck(v, conf-1, vecCa[i], info3);
					fileCheck(v, conf+1, vecCa[i], info4);
					fileCheck(v-1, conf, vecCa[i], info5);
					fileCheck(v+1, conf, vecCa[i], info6);

					if (info1){
						if (i == 0)
							readOutput(n, m, v, conf, vecCa[i], t.data(), si.data());
						else
							readOutput(n, m, v, conf, vecCa[i-1], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
					else if (info2){
						if (i == vecCa.size()-1)
							readOutput(n, m, v, conf, vecCa[i], t.data(), si.data());
						else
							readOutput(n, m, v, conf, vecCa[i+1], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
					else if (info3){
						readOutput(n, m, v, conf-1, vecCa[i], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
					else if (info4){
						readOutput(n, m, v, conf+1, vecCa[i], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
					else if (info5){
						readOutput(n, m, v-1, conf, vecCa[i], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
					else if (info6){
						readOutput(n, m, v+1, conf, vecCa[i], t.data(), si.data());
						for (j = 0; j < m*n; j++){
							sm[j] = si[j];
						}
					}
				}

				// multiple shooting method
				cout << "Shooting for v = " << double(v)/100.0 << ", conf = " 
				     << double(conf)/100.0 << ", Ca = " << Ca << "." << endl;
				mshoot(n, m, nrk, maxiter, tol, Ca, area, vlme, t.data(), si.data(), sf.data(), flag);
			
				// write to file
				if (flag == 0)
					writeSoln(n, m, v, conf, Ca, t.data(), sf.data(), opath);
			
				/* update next initial guess using
				 * first-order continuation */
				if (flag == 0){
					if (i != vecCa.size() - 1){
						if (i == 0) {
							dC0 = vecCa[i  ] - 0;
							dC1 = vecCa[i+1] - 0;
						}
						else {
							dC0 = vecCa[i  ] - vecCa[i-1];
							dC1 = vecCa[i+1] - vecCa[i-1];
						}
						slope = dC1/dC0;
						for (j = 0; j < m*n; j++){
							//si[j] = sm[j] + slope*(sf[j] - sm[j]);
							si[j] = sm[j];
							sm[j] = sf[j];
						}
					}
				}
			}
		}
	}

	return(0);
}



void init(int n, int m, int v, int conf, double Ca,
          double &area, double &vlme, double *t, double *s){
	// declare variables
	int    i, j;
	vector<double> T(m), S(m*n);
	double crit, a;

	// get critical confinement parameter and set nominal radius
	getCritCf(v, crit);
	a = 0.01*double(conf)*crit;
	
	area = 4.0*M_PI*a*a;
	vlme = (4.0/3.0)*M_PI*a*a*a*(0.01*double(v));
	
	// read file
	readEquil(n, m, 90, conf, T.data(), S.data());
	//readOutput(n, m, v, conf-1, Ca, T.data(), S.data()); // initialize using slightly smaller vesicle (useful at high Ca)

	// copy abscissa and solution vectors
	for (i = 0; i < m; i++){
		t[i] = T[i];
		for (j = 0; j < n; j++){
			s[i*n + j] = S[i*n + j];
		}
	}

//	// scale pressure, mean in-plane tension, and transverse shear tension by Ca
//	for (i = 0; i < m; i++){
//		for (j = 3; j < 6; j++){
//			s[i*n + j] = S[i*n + j]/Ca;
//		}
//	}

}



void getCa(bool flip, vector<double>& vecCa){
	// range of capillary numbers
	vecCa.push_back(0.0000001 ); vecCa.push_back(0.00000011); vecCa.push_back(0.00000012); vecCa.push_back(0.00000013); vecCa.push_back(0.00000014);
	vecCa.push_back(0.00000015); vecCa.push_back(0.00000016); vecCa.push_back(0.00000017); vecCa.push_back(0.00000018); vecCa.push_back(0.00000019);
	vecCa.push_back(0.0000002 ); vecCa.push_back(0.00000021); vecCa.push_back(0.00000022); vecCa.push_back(0.00000023); vecCa.push_back(0.00000024);
	vecCa.push_back(0.00000025); vecCa.push_back(0.00000026); vecCa.push_back(0.00000027); vecCa.push_back(0.00000028); vecCa.push_back(0.00000029);
	vecCa.push_back(0.0000003 ); vecCa.push_back(0.00000031); vecCa.push_back(0.00000032); vecCa.push_back(0.00000033); vecCa.push_back(0.00000034);
	vecCa.push_back(0.00000035); vecCa.push_back(0.00000036); vecCa.push_back(0.00000037); vecCa.push_back(0.00000038); vecCa.push_back(0.00000039);
	vecCa.push_back(0.0000004 ); vecCa.push_back(0.00000041); vecCa.push_back(0.00000042); vecCa.push_back(0.00000043); vecCa.push_back(0.00000044);
	vecCa.push_back(0.00000045); vecCa.push_back(0.00000046); vecCa.push_back(0.00000047); vecCa.push_back(0.00000048); vecCa.push_back(0.00000049);
	vecCa.push_back(0.0000005 ); vecCa.push_back(0.00000051); vecCa.push_back(0.00000052); vecCa.push_back(0.00000053); vecCa.push_back(0.00000054);
	vecCa.push_back(0.00000055); vecCa.push_back(0.00000056); vecCa.push_back(0.00000057); vecCa.push_back(0.00000058); vecCa.push_back(0.00000059);
	vecCa.push_back(0.0000006 ); vecCa.push_back(0.00000061); vecCa.push_back(0.00000062); vecCa.push_back(0.00000063); vecCa.push_back(0.00000064);
	vecCa.push_back(0.00000065); vecCa.push_back(0.00000066); vecCa.push_back(0.00000067); vecCa.push_back(0.00000068); vecCa.push_back(0.00000069);
	vecCa.push_back(0.0000007 ); vecCa.push_back(0.00000071); vecCa.push_back(0.00000072); vecCa.push_back(0.00000073); vecCa.push_back(0.00000074);
	vecCa.push_back(0.00000075); vecCa.push_back(0.00000076); vecCa.push_back(0.00000077); vecCa.push_back(0.00000078); vecCa.push_back(0.00000079);
	vecCa.push_back(0.0000008 ); vecCa.push_back(0.00000081); vecCa.push_back(0.00000082); vecCa.push_back(0.00000083); vecCa.push_back(0.00000084);
	vecCa.push_back(0.00000085); vecCa.push_back(0.00000086); vecCa.push_back(0.00000087); vecCa.push_back(0.00000088); vecCa.push_back(0.00000089);
	vecCa.push_back(0.0000009 ); vecCa.push_back(0.00000091); vecCa.push_back(0.00000092); vecCa.push_back(0.00000093); vecCa.push_back(0.00000094);
	vecCa.push_back(0.00000095); vecCa.push_back(0.00000096); vecCa.push_back(0.00000097); vecCa.push_back(0.00000098); vecCa.push_back(0.00000099);

	vecCa.push_back(0.000001 ); vecCa.push_back(0.0000011); vecCa.push_back(0.0000012); vecCa.push_back(0.0000013); vecCa.push_back(0.0000014);
	vecCa.push_back(0.0000015); vecCa.push_back(0.0000016); vecCa.push_back(0.0000017); vecCa.push_back(0.0000018); vecCa.push_back(0.0000019);
	vecCa.push_back(0.000002 ); vecCa.push_back(0.0000021); vecCa.push_back(0.0000022); vecCa.push_back(0.0000023); vecCa.push_back(0.0000024);
	vecCa.push_back(0.0000025); vecCa.push_back(0.0000026); vecCa.push_back(0.0000027); vecCa.push_back(0.0000028); vecCa.push_back(0.0000029);
	vecCa.push_back(0.000003 ); vecCa.push_back(0.0000031); vecCa.push_back(0.0000032); vecCa.push_back(0.0000033); vecCa.push_back(0.0000034);
	vecCa.push_back(0.0000035); vecCa.push_back(0.0000036); vecCa.push_back(0.0000037); vecCa.push_back(0.0000038); vecCa.push_back(0.0000039);
	vecCa.push_back(0.000004 ); vecCa.push_back(0.0000041); vecCa.push_back(0.0000042); vecCa.push_back(0.0000043); vecCa.push_back(0.0000044);
	vecCa.push_back(0.0000045); vecCa.push_back(0.0000046); vecCa.push_back(0.0000047); vecCa.push_back(0.0000048); vecCa.push_back(0.0000049);
	vecCa.push_back(0.000005 ); vecCa.push_back(0.0000051); vecCa.push_back(0.0000052); vecCa.push_back(0.0000053); vecCa.push_back(0.0000054);
	vecCa.push_back(0.0000055); vecCa.push_back(0.0000056); vecCa.push_back(0.0000057); vecCa.push_back(0.0000058); vecCa.push_back(0.0000059);
	vecCa.push_back(0.000006 ); vecCa.push_back(0.0000061); vecCa.push_back(0.0000062); vecCa.push_back(0.0000063); vecCa.push_back(0.0000064);
	vecCa.push_back(0.0000065); vecCa.push_back(0.0000066); vecCa.push_back(0.0000067); vecCa.push_back(0.0000068); vecCa.push_back(0.0000069);
	vecCa.push_back(0.000007 ); vecCa.push_back(0.0000071); vecCa.push_back(0.0000072); vecCa.push_back(0.0000073); vecCa.push_back(0.0000074);
	vecCa.push_back(0.0000075); vecCa.push_back(0.0000076); vecCa.push_back(0.0000077); vecCa.push_back(0.0000078); vecCa.push_back(0.0000079);
	vecCa.push_back(0.000008 ); vecCa.push_back(0.0000081); vecCa.push_back(0.0000082); vecCa.push_back(0.0000083); vecCa.push_back(0.0000084);
	vecCa.push_back(0.0000085); vecCa.push_back(0.0000086); vecCa.push_back(0.0000087); vecCa.push_back(0.0000088); vecCa.push_back(0.0000089);
	vecCa.push_back(0.000009 ); vecCa.push_back(0.0000091); vecCa.push_back(0.0000092); vecCa.push_back(0.0000093); vecCa.push_back(0.0000094);
	vecCa.push_back(0.0000095); vecCa.push_back(0.0000096); vecCa.push_back(0.0000097); vecCa.push_back(0.0000098); vecCa.push_back(0.0000099);

	vecCa.push_back(0.00001 ); vecCa.push_back(0.000011); vecCa.push_back(0.000012); vecCa.push_back(0.000013); vecCa.push_back(0.000014);
	vecCa.push_back(0.000015); vecCa.push_back(0.000016); vecCa.push_back(0.000017); vecCa.push_back(0.000018); vecCa.push_back(0.000019);
	vecCa.push_back(0.00002 ); vecCa.push_back(0.000021); vecCa.push_back(0.000022); vecCa.push_back(0.000023); vecCa.push_back(0.000024);
	vecCa.push_back(0.000025); vecCa.push_back(0.000026); vecCa.push_back(0.000027); vecCa.push_back(0.000028); vecCa.push_back(0.000029);
	vecCa.push_back(0.00003 ); vecCa.push_back(0.000031); vecCa.push_back(0.000032); vecCa.push_back(0.000033); vecCa.push_back(0.000034);
	vecCa.push_back(0.000035); vecCa.push_back(0.000036); vecCa.push_back(0.000037); vecCa.push_back(0.000038); vecCa.push_back(0.000039);
	vecCa.push_back(0.00004 ); vecCa.push_back(0.000041); vecCa.push_back(0.000042); vecCa.push_back(0.000043); vecCa.push_back(0.000044);
	vecCa.push_back(0.000045); vecCa.push_back(0.000046); vecCa.push_back(0.000047); vecCa.push_back(0.000048); vecCa.push_back(0.000049);
	vecCa.push_back(0.00005 ); vecCa.push_back(0.000051); vecCa.push_back(0.000052); vecCa.push_back(0.000053); vecCa.push_back(0.000054);
	vecCa.push_back(0.000055); vecCa.push_back(0.000056); vecCa.push_back(0.000057); vecCa.push_back(0.000058); vecCa.push_back(0.000059);
	vecCa.push_back(0.00006 ); vecCa.push_back(0.000061); vecCa.push_back(0.000062); vecCa.push_back(0.000063); vecCa.push_back(0.000064);
	vecCa.push_back(0.000065); vecCa.push_back(0.000066); vecCa.push_back(0.000067); vecCa.push_back(0.000068); vecCa.push_back(0.000069);
	vecCa.push_back(0.00007 ); vecCa.push_back(0.000071); vecCa.push_back(0.000072); vecCa.push_back(0.000073); vecCa.push_back(0.000074);
	vecCa.push_back(0.000075); vecCa.push_back(0.000076); vecCa.push_back(0.000077); vecCa.push_back(0.000078); vecCa.push_back(0.000079);
	vecCa.push_back(0.00008 ); vecCa.push_back(0.000081); vecCa.push_back(0.000082); vecCa.push_back(0.000083); vecCa.push_back(0.000084);
	vecCa.push_back(0.000085); vecCa.push_back(0.000086); vecCa.push_back(0.000087); vecCa.push_back(0.000088); vecCa.push_back(0.000089);
	vecCa.push_back(0.00009 ); vecCa.push_back(0.000091); vecCa.push_back(0.000092); vecCa.push_back(0.000093); vecCa.push_back(0.000094);
	vecCa.push_back(0.000095); vecCa.push_back(0.000096); vecCa.push_back(0.000097); vecCa.push_back(0.000098); vecCa.push_back(0.000099);

	vecCa.push_back(0.0001 ); vecCa.push_back(0.00011); vecCa.push_back(0.00012); vecCa.push_back(0.00013); vecCa.push_back(0.00014);
	vecCa.push_back(0.00015); vecCa.push_back(0.00016); vecCa.push_back(0.00017); vecCa.push_back(0.00018); vecCa.push_back(0.00019);
	vecCa.push_back(0.0002 ); vecCa.push_back(0.00021); vecCa.push_back(0.00022); vecCa.push_back(0.00023); vecCa.push_back(0.00024);
	vecCa.push_back(0.00025); vecCa.push_back(0.00026); vecCa.push_back(0.00027); vecCa.push_back(0.00028); vecCa.push_back(0.00029);
	vecCa.push_back(0.0003 ); vecCa.push_back(0.00031); vecCa.push_back(0.00032); vecCa.push_back(0.00033); vecCa.push_back(0.00034);
	vecCa.push_back(0.00035); vecCa.push_back(0.00036); vecCa.push_back(0.00037); vecCa.push_back(0.00038); vecCa.push_back(0.00039);
	vecCa.push_back(0.0004 ); vecCa.push_back(0.00041); vecCa.push_back(0.00042); vecCa.push_back(0.00043); vecCa.push_back(0.00044);
	vecCa.push_back(0.00045); vecCa.push_back(0.00046); vecCa.push_back(0.00047); vecCa.push_back(0.00048); vecCa.push_back(0.00049);
	vecCa.push_back(0.0005 ); vecCa.push_back(0.00051); vecCa.push_back(0.00052); vecCa.push_back(0.00053); vecCa.push_back(0.00054);
	vecCa.push_back(0.00055); vecCa.push_back(0.00056); vecCa.push_back(0.00057); vecCa.push_back(0.00058); vecCa.push_back(0.00059);
	vecCa.push_back(0.0006 ); vecCa.push_back(0.00061); vecCa.push_back(0.00062); vecCa.push_back(0.00063); vecCa.push_back(0.00064);
	vecCa.push_back(0.00065); vecCa.push_back(0.00066); vecCa.push_back(0.00067); vecCa.push_back(0.00068); vecCa.push_back(0.00069);
	vecCa.push_back(0.0007 ); vecCa.push_back(0.00071); vecCa.push_back(0.00072); vecCa.push_back(0.00073); vecCa.push_back(0.00074);
	vecCa.push_back(0.00075); vecCa.push_back(0.00076); vecCa.push_back(0.00077); vecCa.push_back(0.00078); vecCa.push_back(0.00079);
	vecCa.push_back(0.0008 ); vecCa.push_back(0.00081); vecCa.push_back(0.00082); vecCa.push_back(0.00083); vecCa.push_back(0.00084);
	vecCa.push_back(0.00085); vecCa.push_back(0.00086); vecCa.push_back(0.00087); vecCa.push_back(0.00088); vecCa.push_back(0.00089);
	vecCa.push_back(0.0009 ); vecCa.push_back(0.00091); vecCa.push_back(0.00092); vecCa.push_back(0.00093); vecCa.push_back(0.00094);
	vecCa.push_back(0.00095); vecCa.push_back(0.00096); vecCa.push_back(0.00097); vecCa.push_back(0.00098); vecCa.push_back(0.00099);

	vecCa.push_back(0.001 ); vecCa.push_back(0.0011); vecCa.push_back(0.0012); vecCa.push_back(0.0013); vecCa.push_back(0.0014);
	vecCa.push_back(0.0015); vecCa.push_back(0.0016); vecCa.push_back(0.0017); vecCa.push_back(0.0018); vecCa.push_back(0.0019);
	vecCa.push_back(0.002 ); vecCa.push_back(0.0021); vecCa.push_back(0.0022); vecCa.push_back(0.0023); vecCa.push_back(0.0024);
	vecCa.push_back(0.0025); vecCa.push_back(0.0026); vecCa.push_back(0.0027); vecCa.push_back(0.0028); vecCa.push_back(0.0029);
	vecCa.push_back(0.003 ); vecCa.push_back(0.0031); vecCa.push_back(0.0032); vecCa.push_back(0.0033); vecCa.push_back(0.0034);
	vecCa.push_back(0.0035); vecCa.push_back(0.0036); vecCa.push_back(0.0037); vecCa.push_back(0.0038); vecCa.push_back(0.0039);
	vecCa.push_back(0.004 ); vecCa.push_back(0.0041); vecCa.push_back(0.0042); vecCa.push_back(0.0043); vecCa.push_back(0.0044);
	vecCa.push_back(0.0045); vecCa.push_back(0.0046); vecCa.push_back(0.0047); vecCa.push_back(0.0048); vecCa.push_back(0.0049);
	vecCa.push_back(0.005 ); vecCa.push_back(0.0051); vecCa.push_back(0.0052); vecCa.push_back(0.0053); vecCa.push_back(0.0054);
	vecCa.push_back(0.0055); vecCa.push_back(0.0056); vecCa.push_back(0.0057); vecCa.push_back(0.0058); vecCa.push_back(0.0059);
	vecCa.push_back(0.006 ); vecCa.push_back(0.0061); vecCa.push_back(0.0062); vecCa.push_back(0.0063); vecCa.push_back(0.0064);
	vecCa.push_back(0.0065); vecCa.push_back(0.0066); vecCa.push_back(0.0067); vecCa.push_back(0.0068); vecCa.push_back(0.0069);
	vecCa.push_back(0.007 ); vecCa.push_back(0.0071); vecCa.push_back(0.0072); vecCa.push_back(0.0073); vecCa.push_back(0.0074);
	vecCa.push_back(0.0075); vecCa.push_back(0.0076); vecCa.push_back(0.0077); vecCa.push_back(0.0078); vecCa.push_back(0.0079);
	vecCa.push_back(0.008 ); vecCa.push_back(0.0081); vecCa.push_back(0.0082); vecCa.push_back(0.0083); vecCa.push_back(0.0084);
	vecCa.push_back(0.0085); vecCa.push_back(0.0086); vecCa.push_back(0.0087); vecCa.push_back(0.0088); vecCa.push_back(0.0089);
	vecCa.push_back(0.009 ); vecCa.push_back(0.0091); vecCa.push_back(0.0092); vecCa.push_back(0.0093); vecCa.push_back(0.0094);
	vecCa.push_back(0.0095); vecCa.push_back(0.0096); vecCa.push_back(0.0097); vecCa.push_back(0.0098); vecCa.push_back(0.0099);

	vecCa.push_back(0.01 ); vecCa.push_back(0.011); vecCa.push_back(0.012); vecCa.push_back(0.013); vecCa.push_back(0.014);
	vecCa.push_back(0.015); vecCa.push_back(0.016); vecCa.push_back(0.017); vecCa.push_back(0.018); vecCa.push_back(0.019);
	vecCa.push_back(0.02 ); vecCa.push_back(0.021); vecCa.push_back(0.022); vecCa.push_back(0.023); vecCa.push_back(0.024);
	vecCa.push_back(0.025); vecCa.push_back(0.026); vecCa.push_back(0.027); vecCa.push_back(0.028); vecCa.push_back(0.029);
	vecCa.push_back(0.03 ); vecCa.push_back(0.031); vecCa.push_back(0.032); vecCa.push_back(0.033); vecCa.push_back(0.034);
	vecCa.push_back(0.035); vecCa.push_back(0.036); vecCa.push_back(0.037); vecCa.push_back(0.038); vecCa.push_back(0.039);
	vecCa.push_back(0.04 ); vecCa.push_back(0.041); vecCa.push_back(0.042); vecCa.push_back(0.043); vecCa.push_back(0.044);
	vecCa.push_back(0.045); vecCa.push_back(0.046); vecCa.push_back(0.047); vecCa.push_back(0.048); vecCa.push_back(0.049);
	vecCa.push_back(0.05 ); vecCa.push_back(0.051); vecCa.push_back(0.052); vecCa.push_back(0.053); vecCa.push_back(0.054);
	vecCa.push_back(0.055); vecCa.push_back(0.056); vecCa.push_back(0.057); vecCa.push_back(0.058); vecCa.push_back(0.059);
	vecCa.push_back(0.06 ); vecCa.push_back(0.061); vecCa.push_back(0.062); vecCa.push_back(0.063); vecCa.push_back(0.064);
	vecCa.push_back(0.065); vecCa.push_back(0.066); vecCa.push_back(0.067); vecCa.push_back(0.068); vecCa.push_back(0.069);
	vecCa.push_back(0.07 ); vecCa.push_back(0.071); vecCa.push_back(0.072); vecCa.push_back(0.073); vecCa.push_back(0.074);
	vecCa.push_back(0.075); vecCa.push_back(0.076); vecCa.push_back(0.077); vecCa.push_back(0.078); vecCa.push_back(0.079);
	vecCa.push_back(0.08 ); vecCa.push_back(0.081); vecCa.push_back(0.082); vecCa.push_back(0.083); vecCa.push_back(0.084);
	vecCa.push_back(0.085); vecCa.push_back(0.086); vecCa.push_back(0.087); vecCa.push_back(0.088); vecCa.push_back(0.089);
	vecCa.push_back(0.09 ); vecCa.push_back(0.091); vecCa.push_back(0.092); vecCa.push_back(0.093); vecCa.push_back(0.094);
	vecCa.push_back(0.095); vecCa.push_back(0.096); vecCa.push_back(0.097); vecCa.push_back(0.098); vecCa.push_back(0.099);

	vecCa.push_back(0.1 ); vecCa.push_back(0.11); vecCa.push_back(0.12); vecCa.push_back(0.13); vecCa.push_back(0.14);
	vecCa.push_back(0.15); vecCa.push_back(0.16); vecCa.push_back(0.17); vecCa.push_back(0.18); vecCa.push_back(0.19);
	vecCa.push_back(0.2 ); vecCa.push_back(0.21); vecCa.push_back(0.22); vecCa.push_back(0.23); vecCa.push_back(0.24);
	vecCa.push_back(0.25); vecCa.push_back(0.26); vecCa.push_back(0.27); vecCa.push_back(0.28); vecCa.push_back(0.29);
	vecCa.push_back(0.3 ); vecCa.push_back(0.31); vecCa.push_back(0.32); vecCa.push_back(0.33); vecCa.push_back(0.34);
	vecCa.push_back(0.35); vecCa.push_back(0.36); vecCa.push_back(0.37); vecCa.push_back(0.38); vecCa.push_back(0.39);
	vecCa.push_back(0.4 ); vecCa.push_back(0.41); vecCa.push_back(0.42); vecCa.push_back(0.43); vecCa.push_back(0.44);
	vecCa.push_back(0.45); vecCa.push_back(0.46); vecCa.push_back(0.47); vecCa.push_back(0.48); vecCa.push_back(0.49);
	vecCa.push_back(0.5 ); vecCa.push_back(0.51); vecCa.push_back(0.52); vecCa.push_back(0.53); vecCa.push_back(0.54);
	vecCa.push_back(0.55); vecCa.push_back(0.56); vecCa.push_back(0.57); vecCa.push_back(0.58); vecCa.push_back(0.59);
	vecCa.push_back(0.6 ); vecCa.push_back(0.61); vecCa.push_back(0.62); vecCa.push_back(0.63); vecCa.push_back(0.64);
	vecCa.push_back(0.65); vecCa.push_back(0.66); vecCa.push_back(0.67); vecCa.push_back(0.68); vecCa.push_back(0.69);
	vecCa.push_back(0.7 ); vecCa.push_back(0.71); vecCa.push_back(0.72); vecCa.push_back(0.73); vecCa.push_back(0.74);
	vecCa.push_back(0.75); vecCa.push_back(0.76); vecCa.push_back(0.77); vecCa.push_back(0.78); vecCa.push_back(0.79);
	vecCa.push_back(0.8 ); vecCa.push_back(0.81); vecCa.push_back(0.82); vecCa.push_back(0.83); vecCa.push_back(0.84);
	vecCa.push_back(0.85); vecCa.push_back(0.86); vecCa.push_back(0.87); vecCa.push_back(0.88); vecCa.push_back(0.89);
	vecCa.push_back(0.9 ); vecCa.push_back(0.91); vecCa.push_back(0.92); vecCa.push_back(0.93); vecCa.push_back(0.94);
	vecCa.push_back(0.95); vecCa.push_back(0.96); vecCa.push_back(0.97); vecCa.push_back(0.98); vecCa.push_back(0.99);

	vecCa.push_back(1.0  );
	vecCa.push_back(2.0  );
	vecCa.push_back(3.0  );
	vecCa.push_back(4.0  );
	vecCa.push_back(5.0  );
	vecCa.push_back(6.0  );
	vecCa.push_back(7.0  );
	vecCa.push_back(8.0  );
	vecCa.push_back(9.0  );

	vecCa.push_back(10.0  );
	vecCa.push_back(11.0  );
	vecCa.push_back(12.0  );
	vecCa.push_back(13.0  );
	vecCa.push_back(14.0  );
	vecCa.push_back(15.0  );
	vecCa.push_back(16.0  );
	vecCa.push_back(17.0  );
	vecCa.push_back(18.0  );
	vecCa.push_back(19.0  );

	vecCa.push_back(20.0  );
	vecCa.push_back(21.0  );
	vecCa.push_back(22.0  );
	vecCa.push_back(23.0  );
	vecCa.push_back(24.0  );
	vecCa.push_back(25.0  );
	vecCa.push_back(26.0  );
	vecCa.push_back(27.0  );
	vecCa.push_back(28.0  );
	vecCa.push_back(29.0  );

	vecCa.push_back(30.0  );
	vecCa.push_back(31.0  );
	vecCa.push_back(32.0  );
	vecCa.push_back(33.0  );
	vecCa.push_back(34.0  );
	vecCa.push_back(35.0  );
	vecCa.push_back(36.0  );
	vecCa.push_back(37.0  );
	vecCa.push_back(38.0  );
	vecCa.push_back(39.0  );

	vecCa.push_back(40.0  );
	vecCa.push_back(41.0  );
	vecCa.push_back(42.0  );
	vecCa.push_back(43.0  );
	vecCa.push_back(44.0  );
	vecCa.push_back(45.0  );
	vecCa.push_back(46.0  );
	vecCa.push_back(47.0  );
	vecCa.push_back(48.0  );
	vecCa.push_back(49.0  );

	vecCa.push_back(50.0  );
	vecCa.push_back(51.0  );
	vecCa.push_back(52.0  );
	vecCa.push_back(53.0  );
	vecCa.push_back(54.0  );
	vecCa.push_back(55.0  );
	vecCa.push_back(56.0  );
	vecCa.push_back(57.0  );
	vecCa.push_back(58.0  );
	vecCa.push_back(59.0  );

	vecCa.push_back(60.0  );
	vecCa.push_back(61.0  );
	vecCa.push_back(62.0  );
	vecCa.push_back(63.0  );
	vecCa.push_back(64.0  );
	vecCa.push_back(65.0  );
	vecCa.push_back(66.0  );
	vecCa.push_back(67.0  );
	vecCa.push_back(68.0  );
	vecCa.push_back(69.0  );

	vecCa.push_back(70.0  );
	vecCa.push_back(71.0  );
	vecCa.push_back(72.0  );
	vecCa.push_back(73.0  );
	vecCa.push_back(74.0  );
	vecCa.push_back(75.0  );
	vecCa.push_back(76.0  );
	vecCa.push_back(77.0  );
	vecCa.push_back(78.0  );
	vecCa.push_back(79.0  );

	vecCa.push_back(80.0  );
	vecCa.push_back(81.0  );
	vecCa.push_back(82.0  );
	vecCa.push_back(83.0  );
	vecCa.push_back(84.0  );
	vecCa.push_back(85.0  );
	vecCa.push_back(86.0  );
	vecCa.push_back(87.0  );
	vecCa.push_back(88.0  );
	vecCa.push_back(89.0  );

	vecCa.push_back(90.0  );
	vecCa.push_back(91.0  );
	vecCa.push_back(92.0  );
	vecCa.push_back(93.0  );
	vecCa.push_back(94.0  );
	vecCa.push_back(95.0  );
	vecCa.push_back(96.0  );
	vecCa.push_back(97.0  );
	vecCa.push_back(98.0  );
	vecCa.push_back(99.0  );

	vecCa.push_back(100.0 );
	vecCa.push_back(110.0 );
	vecCa.push_back(120.0 );
	vecCa.push_back(130.0 );
	vecCa.push_back(140.0 );
	vecCa.push_back(150.0 );
	vecCa.push_back(160.0 );
	vecCa.push_back(170.0 );
	vecCa.push_back(180.0 );
	vecCa.push_back(190.0 );

	vecCa.push_back(200.0 );
	vecCa.push_back(210.0 );
	vecCa.push_back(220.0 );
	vecCa.push_back(230.0 );
	vecCa.push_back(240.0 );
	vecCa.push_back(250.0 );
	vecCa.push_back(260.0 );
	vecCa.push_back(270.0 );
	vecCa.push_back(280.0 );
	vecCa.push_back(290.0 );

	vecCa.push_back(300.0 );
	vecCa.push_back(310.0 );
	vecCa.push_back(320.0 );
	vecCa.push_back(330.0 );
	vecCa.push_back(340.0 );
	vecCa.push_back(350.0 );
	vecCa.push_back(360.0 );
	vecCa.push_back(370.0 );
	vecCa.push_back(380.0 );
	vecCa.push_back(390.0 );

	vecCa.push_back(400.0 );
	vecCa.push_back(410.0 );
	vecCa.push_back(420.0 );
	vecCa.push_back(430.0 );
	vecCa.push_back(440.0 );
	vecCa.push_back(450.0 );
	vecCa.push_back(460.0 );
	vecCa.push_back(470.0 );
	vecCa.push_back(480.0 );
	vecCa.push_back(490.0 );

	vecCa.push_back(500.0 );
	vecCa.push_back(510.0 );
	vecCa.push_back(520.0 );
	vecCa.push_back(530.0 );
	vecCa.push_back(540.0 );
	vecCa.push_back(550.0 );
	vecCa.push_back(560.0 );
	vecCa.push_back(570.0 );
	vecCa.push_back(580.0 );
	vecCa.push_back(590.0 );

	vecCa.push_back(600.0 );
	vecCa.push_back(610.0 );
	vecCa.push_back(620.0 );
	vecCa.push_back(630.0 );
	vecCa.push_back(640.0 );
	vecCa.push_back(650.0 );
	vecCa.push_back(660.0 );
	vecCa.push_back(670.0 );
	vecCa.push_back(680.0 );
	vecCa.push_back(690.0 );

	vecCa.push_back(700.0 );
	vecCa.push_back(710.0 );
	vecCa.push_back(720.0 );
	vecCa.push_back(730.0 );
	vecCa.push_back(740.0 );
	vecCa.push_back(750.0 );
	vecCa.push_back(760.0 );
	vecCa.push_back(770.0 );
	vecCa.push_back(780.0 );
	vecCa.push_back(790.0 );

	vecCa.push_back(800.0 );
	vecCa.push_back(810.0 );
	vecCa.push_back(820.0 );
	vecCa.push_back(830.0 );
	vecCa.push_back(840.0 );
	vecCa.push_back(850.0 );
	vecCa.push_back(860.0 );
	vecCa.push_back(870.0 );
	vecCa.push_back(880.0 );
	vecCa.push_back(890.0 );

	vecCa.push_back(900.0 );
	vecCa.push_back(910.0 );
	vecCa.push_back(920.0 );
	vecCa.push_back(930.0 );
	vecCa.push_back(940.0 );
	vecCa.push_back(950.0 );
	vecCa.push_back(960.0 );
	vecCa.push_back(970.0 );
	vecCa.push_back(980.0 );
	vecCa.push_back(990.0 );

	vecCa.push_back(1000.0 );

	// optional: flip vecCa so largest values appear first
	if (flip){
		double Ca;
		vector<double> vecCaTmp;
		int nCa = vecCa.size();
		int i;
		for (i = 0; i < vecCa.size(); i++){
			Ca = vecCa[nCa-i-1];
			vecCaTmp.push_back(Ca);
		}
		for (i = 0; i < vecCa.size(); i++){
			Ca = vecCaTmp[i];
			vecCa[i] = Ca;
		}
	}
}


void getVr(bool flip, vector<int>& vecVr){
	// range of reduced volumes
//	vecVr.push_back(60);
//	vecVr.push_back(61);
//	vecVr.push_back(62);
//	vecVr.push_back(63);
//	vecVr.push_back(64);
//	vecVr.push_back(65);
//	vecVr.push_back(66);
//	vecVr.push_back(67);
//	vecVr.push_back(68);
//	vecVr.push_back(69);
//
//	vecVr.push_back(70);
	vecVr.push_back(71);
//	vecVr.push_back(72);
//	vecVr.push_back(73);
//	vecVr.push_back(74);
//	vecVr.push_back(75);
//	vecVr.push_back(76);
//	vecVr.push_back(77);
//	vecVr.push_back(78);
//	vecVr.push_back(79);
//	
//	vecVr.push_back(80);
//	vecVr.push_back(81);
//	vecVr.push_back(82);
//	vecVr.push_back(83);
//	vecVr.push_back(84);
//	vecVr.push_back(85);
//	vecVr.push_back(86);
//	vecVr.push_back(87);
//	vecVr.push_back(88);
//	vecVr.push_back(89);
//
//	vecVr.push_back(90);
//	vecVr.push_back(91);
//	vecVr.push_back(92);
//	vecVr.push_back(93);
//	vecVr.push_back(94);
//	vecVr.push_back(95);
//	vecVr.push_back(96);
//	vecVr.push_back(97);
//	vecVr.push_back(98);
//	vecVr.push_back(99);
//	vecVr.push_back(100);

//	// optional: flip vecVr so largest values appear first
//	int nVr = vecVr.size();
//	vector<int> vecVrTmp;
//	for (i = 0; i < vecVr.size(); i++){
//		v = vecVr[nVr-i-1];
//		vecVrTmp.push_back(v);
//	}
//	for (i = 0; i < vecVr.size(); i++){
//		v = vecVrTmp[i];
//		vecVr[i] = v;
//	}

	// optional: flip vecVr so largest values appear first
	if (flip){
		int v, i;
		vector<int> vecVrTmp;
		int nVr = vecVr.size();
		for (i = 0; i < vecVr.size(); i++){
			v = vecVr[nVr-i-1];
			vecVrTmp.push_back(v);
		}
		for (i = 0; i < vecVr.size(); i++){
			v = vecVrTmp[i];
			vecVr[i] = v;
		}
	}
}


void getCf(bool flip, vector<int>& vecCf){
	// range of confinement parameter
	vecCf.push_back(70);
	vecCf.push_back(71);
	vecCf.push_back(72);
//	vecCf.push_back(73);
//	vecCf.push_back(74);
//	vecCf.push_back(75);
//	vecCf.push_back(76);
//	vecCf.push_back(77);
//	vecCf.push_back(78);
//	vecCf.push_back(79);

//	vecCf.push_back(80);
//	vecCf.push_back(81);
//	vecCf.push_back(82);
//	vecCf.push_back(83);
//	vecCf.push_back(84);
//	vecCf.push_back(85);
//	vecCf.push_back(86);
//	vecCf.push_back(87);
//	vecCf.push_back(88);
//	vecCf.push_back(89);
//
//	vecCf.push_back(90);
//	vecCf.push_back(91);
//	vecCf.push_back(92);
//	vecCf.push_back(93);
//	vecCf.push_back(94);
//	vecCf.push_back(95);
//	vecCf.push_back(96);
//	vecCf.push_back(97);
//	vecCf.push_back(98);
//	vecCf.push_back(99);

	// optional: flip vecCf so largest values appear first
	if (flip){
		int nCf = vecCf.size();
		int conf, i;
		vector<int> vecCfTmp;
		for (i = 0; i < vecCf.size(); i++){
			conf = vecCf[nCf-i-1];
			vecCfTmp.push_back(conf);
		}
		for (i = 0; i < vecCf.size(); i++){
			conf = vecCfTmp[i];
			vecCf[i] = conf;
		}
	}

}

void getCritCf(int v, double &crit){
	// set critical confinement parameter
	if (v == 100)
		crit = 1.000000000000000000000000;
	if (v == 99)
		crit = 1.090275107548953613318768;
	if (v == 98)
		crit = 1.133537870150389282270677;
	if (v == 97)
		crit = 1.169545928031930431602916;
	if (v == 96)
		crit = 1.202032009320693326333173;
	if (v == 95)
		crit = 1.232435708492447912665580;
	if (v == 94)
		crit = 1.261494864992458404813790; 
	if (v == 93)
		crit = 1.289649308437067534243365; 
	if (v == 92)
		crit = 1.317187962872965134660298;  
	if (v == 91)
		crit = 1.344314250171293006672972; 
	if (v == 90)
		crit = 1.371179203625788627288043;
	
	if (v == 89)
		crit = 1.397899861949529379050635;
	if (v == 88)
		crit = 1.424570228180523274316834;
	if (v == 87)
		crit = 1.451268169201854679431355;
	if (v == 86)
		crit = 1.478059959100682705404164;
	if (v == 85)
		crit = 1.505003385688966825138370;
	if (v == 84)
		crit = 1.532149944429383541688818;
	if (v == 83)
		crit = 1.559546432719025655752993;
	if (v == 82)
		crit = 1.587236138747886422020625;
	if (v == 81)
		crit = 1.615259749564796639821514;
	if (v == 80)
		crit = 1.643656060707450592693672; 
	
	if (v == 79)
		crit = 1.672462543253015249713266;
	if (v == 78)
		crit = 1.701715807075084193076668;
	if (v == 77)
		crit = 1.731451987830022394460647;
	if (v == 76)
		crit = 1.761707077607608484461251;
	if (v == 75)
		crit = 1.792517213974340291730173; 
	if (v == 74)
		crit = 1.823918938509152593286319;
	if (v == 73)
		crit = 1.855949433369632771788232;
	if (v == 72)
		crit = 1.888646742600595040664463;
	if (v == 71)
		crit = 1.922049983587091034764999;
	if (v == 70)
		crit = 1.956199553113622973856379;  
	
	if (v == 69)
		crit = 1.991137331820519094597204;
	if (v == 68)
		crit = 2.026906890378410524999735;
	if (v == 67)
		crit = 2.063553700384946879653673;
	if (v == 66)
		crit = 2.101125352791323200747707;
	if (v == 65)
		crit = 2.139671786567143914510814;  
	if (v == 64)
		crit = 2.179245530295288117177149;
	if (v == 63)
		crit = 2.219901959443900932366550;
	if (v == 62)
		crit = 2.261699572184745456134028;
	if (v == 61)
		crit = 2.304700286813593219417380;
	if (v == 60)
		crit = 2.348969764079628399293969;
}
