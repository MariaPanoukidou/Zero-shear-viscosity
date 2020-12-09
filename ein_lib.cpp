#include "ein_lib.h"

string in_stress, visc,msd_file;
ifstream infile(in_stress.c_str());
ofstream MSDF(msd_file.c_str());

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

void writeMSDFile(string msd_file, double **msd, double slope1 ,double msdr, int t){
	MSDF.open(msd_file.c_str(), ofstream::out | ofstream::app);
		if(MSDF.is_open()==0){
			cout << "error: could not create output MSD file!" << endl;
			exit(1);
		}
		MSDF.precision(10);
		if(t==0) {MSDF <<  "Ip^2 - Id^2" << setw(20) << "eta" << setw(20) << "msd_r" << setw(20) << "delta t" << endl;}
		MSDF << showpoint  << msd[t][0] << setw(20) << slope1 << setw(20) << msdr  << setw(20) << msd[t][1] << endl;
		MSDF.close();
}

void momentum_msd(double *y,double *x,double **msd, int Nconfs){
	int numberOfdeltaT = Nconfs;

		 for(int dt = 1; dt<numberOfdeltaT;dt++){
			 vector<double> deltaCoords(numberOfdeltaT,0.0);
			 for(int i = 0; i<numberOfdeltaT-dt-1;i++){
				 
					deltaCoords[dt-1] += (y[i+dt] - y[i])*(x[i+dt] - x[i]);
					
			 }
			 
			 deltaCoords[dt-1] /= (double) (numberOfdeltaT-dt);

           msd[dt-1][0] =  deltaCoords[dt-1] ;
           msd[dt-1][1] = dt;
		  
		 }
}

void trapz(double *I1xy,double *I1xz,double *I1yz, double *sxy1,double *sxz1, double *syz1,double *I2xy,double *I2xz, double *I2yz, double *sxy2,double *sxz2, double *syz2,double *I3xy,double *I3xz,double *I3yz, double *sxy3,double *sxz3, double *syz3,double *Time, int N){
			// trapezoidal rule over the stress values
			int j = 1;
			double dt,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9;
			while (j - 1 <= N - 1) {
				sum1 = 0.0;sum2 = 0.0;sum3 = 0.0;
				sum4 = 0.0;sum5 = 0.0;sum6 = 0.0;
				sum7 = 0.0;sum8 = 0.0;sum9 = 0.0;
				for (int i = 0; i < j; i++) {
					dt = Time[i+1]-Time[i];
					sum1 += dt*((sxy1[i+1]+sxy1[i])/2.0);
					sum2 += dt*((sxz1[i+1]+sxz1[i])/2.0);
					sum3 += dt*((syz1[i+1]+syz1[i])/2.0);
					
					sum4 += dt*((sxy2[i+1]+sxy2[i])/2.0);
					sum5 += dt*((sxz2[i+1]+sxz2[i])/2.0);
					sum6 += dt*((syz2[i+1]+syz2[i])/2.0);
					
					sum7 += dt*((sxy3[i+1]+sxy3[i])/2.0);
					sum8 += dt*((sxz3[i+1]+sxz3[i])/2.0);
					sum9 += dt*((syz3[i+1]+syz3[i])/2.0);
					
				}
				I1xy[j] = sum1;
				I1xz[j] = sum2;
				I1yz[j] = sum3;
				
				I2xy[j] = sum4;
				I2xz[j] = sum5;
				I2yz[j] = sum6;

				I3xy[j] = sum7;
				I3xz[j] = sum8;
				I3yz[j] = sum9;
				j++;
			}	
}