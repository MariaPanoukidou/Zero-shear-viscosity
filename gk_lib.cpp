#include "gk_lib.h"

string in_stress, in_stress2, stress_d, stress_k, stress_r, AC, visc;
ifstream infile(in_stress.c_str());
ifstream infile_2(in_stress2.c_str());
ifstream infile_d(stress_d.c_str());
ifstream infile_r(stress_r.c_str());
ofstream autocorrF(AC.c_str());
ofstream viscF(visc.c_str());

void ReadStressFile(string in_stress, double *sxy,double *sxz, double *syz, double *t, int Nconfs){

	istringstream iss;
	string line;
	infile.open(in_stress.c_str());
	 if (!infile)
    	{
        // Print an error and exit
        cerr << "input Stress file could not be opened for reading!" << endl;
        exit(1);
    	}
	 for(int i=0;i<1;i++){
		getline(infile,line);
	 }
	int Line = 0;
	double tr, t0, t1, press, sxx, syx, syy, szx, szy, szz, temp;
	while(!infile.eof()){
			getline(infile,line);
			iss.clear();
			iss.str(line);
			if(Line == 0) iss >> t0;
			else if(Line == 1) iss >> t1;
			else iss >> tr;
			iss >> press;
			iss >> sxx;
			iss >> sxy[Line];
			iss >> sxz[Line];
			iss >> syx;
			iss >> syy;
			iss >> syz[Line];
			iss >> szx;
			iss >> szy;
			iss >> szz;
			iss >> temp;
			Line++;
		}
	infile.close();

	for(int k = 0; k < Nconfs; k++){
		t[k] = k*(t1-t0) + (t1-t0);
	}
}

void writeAutocorrFile(string AC, double *mtoxy,double *mtoxz, double *mtoyz, double *tt, int t){
	autocorrF.open(AC.c_str(), ofstream::out | ofstream::app);
		if(autocorrF.is_open()==0){
			cout << "error: could not create output autocorrelation file!" << endl;
			exit(1);
		}
		autocorrF.precision(10);
		if(t==0){
			autocorrF << "time" << setw(20) << "ACxy" << setw(20) << "ACxz" << setw(20) << "ACyz" << endl;
			autocorrF << showpoint << tt[t] << setw(20) << stoxy[t] << setw(20) << stoxz[t] << setw(20) << stoyz[t] << endl;
			autocorrF.close();
		}else{
			autocorrF << showpoint << tt[t] << setw(20) << stoxy[t] << setw(20) << stoxz[t] << setw(20) << stoyz[t] << endl;
			autocorrF.close();
		}
}

void writeViscFile(string visc, double etaxy,double etaxz,double etayz, double *tt, int t){
	viscF.open(visc.c_str(), ofstream::out | ofstream::app);
		if(viscF.is_open()==0){
			cout << "error: could not create output viscosity file!" << endl;
			exit(1);
		}
		viscF.precision(10);
		if(t==0) viscF << "time" << setw(20) << "viscosity_xy" << setw(20) << "viscosity_xz" << setw(20) << "viscosity_yz" << endl;
		viscF << showpoint << tt[t] << setw(20) << etaxy << setw(20) << etaxz << setw(20) << etayz << endl;
		viscF.close();

}

double CalcAutocorr(double *mtoValue, double *s, double *sd, double *sr, int i, int Nconfs){
		double AcR0;		
		int nmax = Nconfs - i;
		double Valuec;
		double Valued;		
		double Valuecd;
		double Valuedc;
		double ValueR;
		Valuec = 0.0;
		Valued = 0.0;
		Valuecd = 0.0;
		Valuedc = 0.0;
		ValueR = 0.0;
		for(int z = 0; z<nmax; z++){
			Valuec = Valuec + s[i+z]*s[z];
			Valued = Valued + sd[i+z]*sd[z];
			Valuecd = Valuecd + sd[i+z]*s[z];
			Valuedc = Valuedc + s[i+z]*sd[z];
			ValueR = ValueR + sr[i+z]*sr[z];
		}
		mtoValue[i] = Valuec + Valuecd - Valuedc - Valued;
		mtoValue[i] = mtoValue[i]/(double) Nconfs;
		AcR0 = ValueR/(double) Nconfs;
	return AcR0;
}
