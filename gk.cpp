//This code was developed by M. Panoukidou (University of Manchester) ( February 2020 )
// within the frames of her PhD research and for the purposes of the work:
//
//              M. Panoukidou, CR Wand, P. Carbobe, Phys.Rew E (2020) 
//
// The Stress Autocorrelation function is calculated using multiple time origins (MTO).
// Afterwards the viscosity is calculated using the generalized Green- Kubo relation for DPD derived by: 
//
// 				M.H. Ernst, R. Brito, Phys. Rev. E. 72 (2005) 1–11
// 				M.H. Ernst, R. Brito, Europhys. Lett. 73 (2006) 183
// 				G. Jung, F. Schmid, J. Chem. Phys. 144 (2016) 204104
//
// If you are using this code please acknowledge our effort and cite our paper "M. Panoukidou, CR Wand, P. Carbobe, JCTC (2020)".
// Note: This cpp file needs its headers to run. Please compile with the header file gk_lib.h and gk_lib.cpp.

#include "gk_lib.h"

using namespace std;
int Nconfs, N;
double V;
double kbT;

int main(int argc, char *argv[])
{
	if (argc!=5){
		cout << "error: wrong invocation!" << endl;
		cout << "try with:" << endl;
		cout << " distribution <int Nconfs> <double Volume> <double kbT> <int Nintegral>" << endl;
		return 0; 
	}
	
	Nconfs = atoi(argv[1]);
	V = atof(argv[2]);
	kbT = atof(argv[3]);
	N = atoi(argv[4]);
	
	string input_file1,input_file2,input_file3,input_file4;
	ifstream infile;	
	ofstream outfile,outfile2,outfile3;

	double *sxy = NULL;
	double *sxz = NULL;
	double *syz = NULL;
	
	double *scxy = NULL;
	double *scxz = NULL;
	double *scyz = NULL;
	
	double *srxy = NULL;
	double *srxz = NULL;
	double *sryz = NULL;
	
	double *sdxy = NULL;
	double *sdxz = NULL;
	double *sdyz = NULL;
	
	double *skxy = NULL;
	double *skxz = NULL;
	double *skyz = NULL;
	
	double *MTOxy = NULL;
	double *MTOxz = NULL;
	double *MTOyz = NULL;
	
	double *Time = NULL;

	input_file1 = "Stress_pot.d";
	input_file2 = "Stress_diss.d";
	input_file3 = "Stress_rn.d";
	input_file4 = "Stress_kin.d";
	AC = "pressure_autocorr.dat";
	visc = "viscosity_integral.dat";

	double AcR0_xy,AcR0_xz,AcR0_yz;

	sxy = new double[Nconfs];
	sxz = new double[Nconfs];
	syz = new double[Nconfs];
	
	scxy = new double[Nconfs];
	scxz = new double[Nconfs];
	scyz = new double[Nconfs];
	
	srxy = new double[Nconfs];
	srxz = new double[Nconfs];
	sryz = new double[Nconfs];
	
	sdxy = new double[Nconfs];
	sdxz = new double[Nconfs];
	sdyz = new double[Nconfs];
	
	skxy = new double[Nconfs];
	skxz = new double[Nconfs];
	skyz = new double[Nconfs];
	
	Time = new double[Nconfs];
	MTOxy = new double[Nconfs];
	MTOxz = new double[Nconfs];
	MTOyz = new double[Nconfs];

// initialization of arrays 
	for(int i = 0; i < Nconfs; i++){
		
		MTOxy[i] = 0.0;
		MTOxz[i] = 0.0;
		MTOyz[i] = 0.0;
		
		sxy[i] = 0.0;
		sxz[i] = 0.0;
		syz[i] = 0.0;
		
		scxy[i] = 0.0;
		scxz[i] = 0.0;
		scyz[i] = 0.0;
		
		sdxy[i] = 0.0;
		sdxz[i] = 0.0;
		sdyz[i] = 0.0;
		
		srxy[i] = 0.0;
		srxz[i] = 0.0;
		sryz[i] = 0.0;
		
		skxy[i] = 0.0;
		skxz[i] = 0.0;
		skyz[i] = 0.0;
		
		Time[i] = 0.0;
	}

	ReadStressFile(input_file1,scxy,scxz,scyz,Time,Nconfs);
	ReadStressFile(input_file2,sdxy,sdxz,sdyz,Time,Nconfs);
	ReadStressFile(input_file3,srxy,srxz,sryz,Time,Nconfs);
	ReadStressFile(input_file4,skxy,skxz,skyz,Time,Nconfs);
	cout << "I have read all stress files!" << endl;
	for(int j = 0; j<Nconfs; j++){
		sxy[j] = scxy[j] + skxy[j];
		sxz[j] = scxz[j] + skxz[j];
		syz[j] = scyz[j] + skyz[j];
	}

	double r0xy, r0xz, r0yz;
		//Calculate STO
	for(int i = 0; i<N; i++){
		r0xy = CalcAutocorr(MTOxy,sxy,sdxy,srxy,i,N);
		r0xz = CalcAutocorr(MTOxz,sxz,sdxz,srxz,i,N);
		r0yz = CalcAutocorr(MTOyz,syz,sdyz,sryz,i,N);
		if(i==0){
			AcR0_xy = r0xy;
			AcR0_xz = r0xz;
			AcR0_yz = r0yz;
		} 

		//write STO into a file in order to plot	
		writeAutocorrFile(AC,MTOxy,MTOxz,MTOyz,Time,i);
	}
	

//Calculate the viscosity using the Green - Kubo relation (integrate by trapezoidal algorithm)
		vector<double> Ixy(N,0.0);
		vector<double> Ixz(N,0.0);
		vector<double> Iyz(N,0.0);
  		int j = 1;
  		double sum1, sum2, sum3, etaxy, etaxz, etayz, dt;
  		while (j - 1 <= N - 1) {
    		sum1 = 0.0;sum2 = 0.0; sum3 = 0.0;
    		for (int i = 0; i < j; i++) {
      			dt = Time[i+1]-Time[i];
      			sum1 += dt*((MTOxy[i+1]+MTOxy[i])/2.0);
				sum2 += dt*((MTOxz[i+1]+MTOxz[i])/2.0);
				sum3 += dt*((MTOyz[i+1]+MTOyz[i])/2.0);
    		}

			Ixy[j] = sum1 + (dt/2.0)*AcR0_xy;
			Ixz[j] = sum2 + (dt/2.0)*AcR0_xz;
			Iyz[j] = sum3 + (dt/2.0)*AcR0_yz;

			Ixy[j] *= V;
			Ixz[j] *= V;
			Iyz[j] *= V;
			Ixy[j] /= kbT;
			Ixz[j] /= kbT;
			Iyz[j] /= kbT;

			writeViscFile(visc,Ixy[j],Ixz[j],Iyz[j],Time,j-1);	
    		j++;
  		}
	
	delete [] MTOxy;
	delete [] MTOxz;
	delete [] MTOyz;
 	delete [] Time;
	delete [] sxy;
	delete [] sxz;
	delete [] syz;
	delete [] scxy;
	delete [] scxz;
	delete [] scyz;
	delete [] sdxy;
	delete [] sdxz;
	delete [] sdyz;
	delete [] srxy;
	delete [] srxz;
	delete [] sryz;
	delete [] skxy;
	delete [] skxz;
	delete [] skyz;
	return 0;
}

