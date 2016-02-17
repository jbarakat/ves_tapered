/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
//#include <algorithm>
//#include <iterator>
#include <string>
#include <gsl/gsl_sf_trig.h>
#include "../include/shoot.h"
#include "../include/read.h"
#include "../include/write.h"

void init(int, int, int, int, double &, double&, double&, double*, double*);
void getCa(bool, vector<double>&);
void getVr(bool, vector<int   >&);
void getTa(bool, vector<int   >&);
void getCf(bool, vector<int   >&);
void getCritCf(int, double&);

int main(){
	int    i, j, k, l;
	int    n    = 12;		// size of solution vector
	int    m    = 105;	// number of shooting points
	int    nrk  = 50;		// number of Runge-Kutta steps
	int    flag;				// error flag
	bool   info;				// boolean flag for file check

	// output directory
	string opath = "../output";

	// abscissa and solution vectors
	vector<double> t(m), si(n*m), sm(n*m), sf(n*m);

	// parameters for time integration
	double ts, dts;
	double r0, r1, dr;
	double x0, x1, dx;
	double cos, sin;
	double U;
	vector<double> u(m, 0.0);
	
	// parameters
	double par[10];
	double Ca, area, vlme, R0, alph, tana;
	double redvol, tubrad, tapang;
	double dC0, dC1, slope;
	vector<double> vecCa;
	vector<int   > vecVr;
	vector<int   > vecCf;
	vector<int   > vecTa;
	
	int    v;						// reduced volume
	int    conf;				// initial confinement
	int    a;						// taper angle
	int    id[10];			// identifier for file saving

	// initialize Ca
	Ca = 10000;

	// get vector of capillary numbers
	bool flipCa = false;
	getCa(flipCa, vecCa);

	// get vector of reduced volumes
	bool flipVr = true;
	getVr(flipVr, vecVr);
	
	// get vector of taper angle
	bool flipTa = false;
	getTa(flipTa, vecTa);
	tapang = 0.01;

	// get vector of confinement
	bool flipCf = false;
	getCf(flipCf, vecCf);

	// loop over reduced volume
	for (l = 0; l < vecVr.size(); l++){
		v = vecVr[l];

		// loop over confinement
		for (k = 0; k < vecCf.size(); k++){
			conf = vecCf[k];

			// loop over capillary number to check
			// which input file to use
			for (i = 0; i < vecCa.size(); i++){
				Ca = vecCa[i];
				if (conf < 100)
					fileCheckInput(v, conf, Ca, info);
				else
					fileCheckInput(v, 99, Ca, info);

				if (info){ // file exists
					break;
				}
				Ca = 0.0;
			}
			
			if (fabs(Ca) < 1e-12){
				cout << "No input file... Skipping to next v, conf..." << endl;
				continue;
			}

			// initialize
			if (l == 0 && k == 0){
				cout << "Initializing... " << endl;
			}
			init(n, m, v, conf, Ca, redvol, tubrad, t.data(), si.data());

			/* update next initial guess using
			 * numerical continuation */
			if (l > 0 || k > 0){	
				if (flag == 0){
					cout << "Continuing from last solution..." << endl;
					for (j = 0; j < m*n; j++){
						si[j] = sf[j];
					}
				}
			}

			// set parameters
			par[0] = redvol;
			par[1] = tubrad;
			par[2] = tapang;

			/*--------------------- TIME INTEGRATION ----------------------*/
		
			
			ts   = 0.0 ;
			dts  = 10;
			R0   = tubrad;
			alph = tapang;
			tana = gsl_sf_sin(alph)/gsl_sf_cos(alph);
			
			// initialize to get the first value of u
			vector<double> s0;
			for (i = 0; i < 3; i++){
				cout << "Shooting for v = " << 0.01*v << ", conf = 0."  << conf << "." << endl;
				mshoot(n, m, nrk, par, u.data(), t.data(), si.data(), sf.data(), flag);

				if (flag == 0){
					cout << "Continuing from last solution..." << endl;
					for (j = 0; j < m*n; j++){
						si[j] = sf[j];
					}
				}
			}

			// solve quasi-steady problem using the multiple shooting method
			for (i = 0; i < 100; i++){
				cout << "Shooting for v = " << 0.01*v << ", conf = 0."  << conf << "." << endl;
				mshoot(n, m, nrk, par, u.data(), t.data(), si.data(), sf.data(), flag);
				
				// update R0
				U      = sf[0*n + 8];
				R0    += -U*dts*tana;
				par[1] = R0;

				// START FROM HERE

				// calculate source field for next timestep
				// (normal velocity of the surface)
				if (i > 4){
					double unorm = 0.0;
					for (j = 0; j < m; j++){
						cos  = gsl_sf_cos(si[j*n + 1]);
						sin  = gsl_sf_sin(si[j*n + 1]);
						r0   = si[j*n + 0 ];
						r1   = sf[j*n + 0 ];
						x0   = si[j*n + 10];
						x1   = sf[j*n + 10];
						dr   = (r1 - r0)/dts;
						dx   = (x1 - x0)/dts;

						u[j] = dx*sin + dr*cos;
						unorm += u[j]*u[j];
					}
					cout << unorm << endl;
				}
			
				// write to file
				int a = 1;
				id[0] = v;
				id[1] = 1;
				if (flag == 0)
					writeSoln(n, m, id, t.data(), sf.data(), opath);

				if (flag == 0){
					cout << "Continuing from last solution..." << endl;
					for (j = 0; j < m*n; j++){
						si[j] = sf[j];
					}
				}

				ts += dts;
			}
			
			/*-------------------------------------------------------------*/
		}
	}

	return(0);
}






void init(int n, int m, int v, int conf, double &Ca,
          double &redvol, double &tubrad, double *t, double *s){
	// declare variables
	int    i, j;
	double crit, area, vlme;
	vector<double> T(m), S(m*n);

	// get critical confinement for given reduced volume
	getCritCf(v, crit);

	// set nominal radius
	double a, a2, a3;
	if (conf < 100)
		a = 0.01*double(conf)*crit;
	else 
		a = 0.001*double(conf)*crit;
	
	// set area and volume
	area = 4.0*M_PI*a*a;
	vlme = (4.0/3.0)*M_PI*a*a*a*(0.01*double(v));
	
	// read file
	if (conf < 100)
		readInput(n, m, v, conf, Ca, T.data(), S.data());
	else
		readInput(n, m, v, 99, Ca, T.data(), S.data());

	// copy abscissa and solution vectors
	for (i = 0; i < m; i++){
		t[i] = T[i];
		for (j = 0; j < n; j++){
			s[i*n + j] = S[i*n + j];
		}
	}

	/* NOTE: The variables read in from the input files are scaled by
	 * the bending modulus, tube radius, and suspending fluid viscosity.
	 * Need to switch to a set of dimensionless variables that use the
	 * pressure drop, nominal vesicle radius, and suspending fluid viscosity
	 * as the characteristic scales (see below). */

	// scale pressure and surface tension by Ca (this uses
	// the linear velocity instead of the bending modulus
	// as the characteristic scale)
	for (i = 0; i < m; i++){
		double p   = s[i*n + 2]/Ca;
		double sig = s[i*n + 3]/Ca;

		s[i*n + 2] = p  ;
		s[i*n + 3] = sig;
	}
	
	// rescale all variables by pressure drop
	double dp = s[0 + 2] - s[(m-1)*n + 2];
	for (i = 0; i < m; i++){
		s[i*n + 2 ] /= dp; // p
		s[i*n + 3 ] /= dp; // sig
		s[i*n + 6 ] /= dp; // Q2
		s[i*n + 8 ] /= dp; // U
	}

	// rescale all variables by the nominal radius of the vesicle
	a  = sqrt(area/(4.0*M_PI));
	a2 = a*a;
	a3 = a*a*a;
	for (i = 0; i < m; i++){
		s[i*n + 0 ] /= a ; // r
		s[i*n + 2 ] *= a ; // p
		s[i*n + 4 ] /= a2; // A 
		s[i*n + 5 ] /= a3; // V
		s[i*n + 6 ] /= a ; // Q2
		s[i*n + 7 ] /= a ; // Rt
		s[i*n + 9 ] /= a ; // S
		s[i*n + 10] /= a ; // x
		s[i*n + 11] /= a ; // xcm
	}
	
	// get tube radius and reduced volume
	redvol = vlme/(4.0*M_PI*a3/3.0);
	tubrad = 1.0/a;

	// scale area and volume wrt nominal radius of vesicle
	area = area/a2;
	vlme = vlme/a3;
}









void getCa(bool flip, vector<double>& vecCa){
	// range of capillary number (for input file only)
	vecCa.push_back(1000.0   );
	vecCa.push_back( 100.0   );
	vecCa.push_back(  10.0   );
	vecCa.push_back(   1.0   );
	vecCa.push_back(   0.0010);
}

void getTa(bool flip, vector<int>& vecTa){
	// range of taper angle
	vecTa.push_back(0);
	vecTa.push_back(1);
	vecTa.push_back(2);
	vecTa.push_back(3);
	vecTa.push_back(4);
	vecTa.push_back(5);
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

//	vecVr.push_back(70);
//	vecVr.push_back(71);
//	vecVr.push_back(72);
//	vecVr.push_back(73);
//	vecVr.push_back(74);
//	vecVr.push_back(75);
//	vecVr.push_back(76);
//	vecVr.push_back(77);
//	vecVr.push_back(78);
//	vecVr.push_back(79);
	
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

	vecVr.push_back(90);
//	vecVr.push_back(91);
//	vecVr.push_back(92);
//	vecVr.push_back(93);
//	vecVr.push_back(94);
//	vecVr.push_back(95);
//	vecVr.push_back(96);
//	vecVr.push_back(97);
//	vecVr.push_back(98);
//	vecVr.push_back(99);
//	//vecVr.push_back(100);

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
//	vecCf.push_back(70);
//	vecCf.push_back(71);
//	vecCf.push_back(72);
//	vecCf.push_back(73);
//	vecCf.push_back(74);
//	vecCf.push_back(75);
//	vecCf.push_back(76);
//	vecCf.push_back(77);
//	vecCf.push_back(78);
//	vecCf.push_back(79);

	vecCf.push_back(80);
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
//
//	vecCf.push_back(991);
//	vecCf.push_back(992);
//	vecCf.push_back(993);
//	vecCf.push_back(994);
//	vecCf.push_back(995);
//	vecCf.push_back(996);
//	vecCf.push_back(997);
//	vecCf.push_back(998);
//	vecCf.push_back(999);

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
