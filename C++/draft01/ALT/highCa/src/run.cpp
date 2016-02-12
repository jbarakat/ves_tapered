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
#include "../include/shoot.h"
#include "../include/read.h"
#include "../include/write.h"

void init(int, int, int, int, double, double&, double&, double*, double*);

int main(){
	int    i, j, k, l;
	int    n    = 10;		// size of solution vector
	int    m    = 105;	// number of shooting points
	int    nrk;					// number of Runge-Kutta steps
	int    v;						// reduced volume
	int    conf;				// maximum radius
	int    flag;				// error flag
	bool   info;				// boolean flag for file check

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

	// initialize Ca
	Ca = 10000;

	// range of capillary number (for input file only)
	vecCa.push_back(1000.0   );
	vecCa.push_back( 100.0   );
	vecCa.push_back(  10.0   );
	vecCa.push_back(   1.0   );
	vecCa.push_back(   0.0010);

	// range of reduced volumes
//	vecVr.push_back(60);
//	vecVr.push_back(61);
//	vecVr.push_back(62);
//	vecVr.push_back(63);
//	vecVr.push_back(64);
	vecVr.push_back(65);
//	vecVr.push_back(66);
//	vecVr.push_back(67);
//	vecVr.push_back(68);
//	vecVr.push_back(69);

	vecVr.push_back(70);
//	vecVr.push_back(71);
//	vecVr.push_back(72);
//	vecVr.push_back(73);
//	vecVr.push_back(74);
	vecVr.push_back(75);
//	vecVr.push_back(76);
//	vecVr.push_back(77);
//	vecVr.push_back(78);
//	vecVr.push_back(79);
	
	vecVr.push_back(80);
//	vecVr.push_back(81);
//	vecVr.push_back(82);
//	vecVr.push_back(83);
//	vecVr.push_back(84);
	vecVr.push_back(85);
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
	int nVr = vecVr.size();
	vector<int> vecVrTmp;
	for (i = 0; i < vecVr.size(); i++){
		v = vecVr[nVr-i-1];
		vecVrTmp.push_back(v);
	}
	for (i = 0; i < vecVr.size(); i++){
		v = vecVrTmp[i];
		vecVr[i] = v;
	}

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
	vecCf.push_back(81);
	vecCf.push_back(82);
	vecCf.push_back(83);
	vecCf.push_back(84);
	vecCf.push_back(85);
	vecCf.push_back(86);
	vecCf.push_back(87);
	vecCf.push_back(88);
	vecCf.push_back(89);

	vecCf.push_back(90);
	vecCf.push_back(91);
	vecCf.push_back(92);
	vecCf.push_back(93);
	vecCf.push_back(94);
	vecCf.push_back(95);
	vecCf.push_back(96);
	vecCf.push_back(97);
	vecCf.push_back(98);
	vecCf.push_back(99);

	vecCf.push_back(991);
	vecCf.push_back(992);
	vecCf.push_back(993);
	vecCf.push_back(994);
	vecCf.push_back(995);
	vecCf.push_back(996);
	vecCf.push_back(997);
	vecCf.push_back(998);
	vecCf.push_back(999);

//	// optional: flip vecCf so largest values appear first
//	int nCf = vecCf.size();
//	vector<int> vecCfTmp;
//	for (i = 0; i < vecCf.size(); i++){
//		conf = vecCf[nCf-i-1];
//		vecCfTmp.push_back(conf);
//	}
//	for (i = 0; i < vecCf.size(); i++){
//		conf = vecCfTmp[i];
//		vecCf[i] = conf;
//	}
	

	
	// loop over reduced volume
	for (l = 0; l < vecVr.size(); l++){
		v = vecVr[l];

		// loop over confinement parameter
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
			nrk = 50;
			cout << "Initializing... " << endl;
			init(n, m, v, conf, Ca, area, vlme, t.data(), si.data());
			for (j = 0; j < m*n; j++){
				sm[j] = si[j];
			}
			cout << "Initialization complete." << endl;
			
//			// check neighboring solutions
//			fileCheckOutput(v, (conf-1)/10, info); // for very high confinement, e.g. conf = 991
//			if (info){ // file exists
//				readOutput(n, m, v, (conf-1)/10, t.data(), si.data());
//				for (j = 0; j < m*n; j++){
//					sm[j] = si[j];
//				}
//			}
//
//			fileCheckOutput(v, conf-1, info);
//			if (info){ // file exists
//				readOutput(n, m, v, conf-1, t.data(), si.data());
//				for (j = 0; j < m*n; j++){
//					sm[j] = si[j];
//				}
//			}



			
		//	fileCheckOutput(v, conf+1, info);
		//	if (info){ // file exists
		//		readOutput(n, m, v, conf+1, t.data(), si.data());
		//		for (j = 0; j < m*n; j++){
		//			sm[j] = si[j];
		//		}
		//	}
			
		//	fileCheckOutput(v-1, conf, info);
		//	if (info){ // file exists
		//		readOutput(n, m, v-1, conf, t.data(), si.data());
		//		for (j = 0; j < m*n; j++){
		//			sm[j] = si[j];
		//		}
		//	}
			
		//	fileCheckOutput(v+1, conf, info);
		//	if (info){ // file exists
		//		readOutput(n, m, v+1, conf, t.data(), si.data());
		//		for (j = 0; j < m*n; j++){
		//			sm[j] = si[j];
		//		}
		//	}

		//	fileCheckOutput(v, conf, info); 
		//	if (info){ // file exists
		//		readOutput(n, m, v, conf, t.data(), si.data());
		//		for (j = 0; j < m*n; j++){
		//			sm[j] = si[j];
		//		}
		//	}

			// multiple shooting method
			cout << "Shooting for v = " << 0.01*v << ", conf = 0." 
			     << conf << ", Ca = " << Ca << "." << endl;
			mshoot(n, m, nrk, Ca, area, vlme, t.data(), si.data(), sf.data(), flag);
			
//			// write to file
//			if (flag == 0)
//				writeSoln(n, m, v, conf, t.data(), sf.data(), opath);
		}
	}

	return(0);
}

void init(int n, int m, int v, int conf, double Ca,
          double &area, double &vlme, double *t, double *s){
	// declare variables
	int    i, j;
	double crit;
	vector<double> T(m), S(m*n);

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

	// set nominal radius
	double a;
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

	// scale pressure and surface tension by Ca
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

//	for (i = 0; i < m; i++){
//		for (j = 0; j < n; j++){
//		}
//	}
}

