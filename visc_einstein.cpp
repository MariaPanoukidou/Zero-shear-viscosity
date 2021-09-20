//This code was developed by M. Panoukidou (University of Manchester) ( February 2020 )
//The zero-shear viscosity is calculated using the generalized Eintein-Helfand formula derived in the work:
//
//              M. Panoukidou, CR Wand, P. Carbobe, Soft Matter (2021), DOI: 10.1039/d1sm00891a
//
// If you are using this code please acknowledge our effort and cite our paper "M. Panoukidou, CR Wand, P. Carbobe, DOI: 10.1039/d1sm00891a".
// Note: This cpp file needs its headers to run. Please compile with the header file ein_lib.h and ein_lib.cpp.

#include "ein_lib.h"

using namespace std;
int Nconfs, N, rStart;
double V;
double kbT;

int main(int argc, char *argv[])
{
	if (argc!=5){
		cout << "error: wrong invocation!" << endl;
		cout << "try with:" << endl;
		cout << " distribution <int Nconfs> <double Volume> <double kbT> <int Nintegral> " << endl;
		return 0; 
	}
	
	Nconfs = atoi(argv[1]);
	V = atof(argv[2]);
	kbT = atof(argv[3]);
	N = atoi(argv[4]);
	
	string input_file1,input_file2,input_file3,input_file4;
	ifstream infile,readstoF;	
	ofstream outfile,outfile2,outfile3, restart_file;

	double *spxz = NULL;
	double *spxy = NULL;
	double *spyz= NULL;
	
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
	
	double *Ipxy = NULL;
	double *Ipxz = NULL;
	double *Ipyz = NULL;
	double *Idxy = NULL;
	double *Idxz = NULL;
	double *Idyz = NULL;
	double *Irxy  = NULL;
	double *Irxz  = NULL;
	double *Iryz  = NULL;

	double *Time = NULL;
	
	double **msd;
	double **msdp;
	double **msdd;
	double **msdr;
	
	double **msd_pxy;
	double **msd_pxz;
	double **msd_pyz;

	double **msd_dxy;
	double **msd_dxz;
	double **msd_dyz;
	
	double **msd_rxy;
	double **msd_rxz;
	double **msd_ryz;

	input_file1 = "Stress_pot.d";
	input_file2 = "Stress_diss.d";
	input_file3 = "Stress_rn.d";
	input_file4 = "Stress_kin.d";
	msd_file = "momentum_diff.dat";

	double AcR0_xy,AcR0_xz,AcR0_yz;

	// memory allocation
	spxz = new double[Nconfs];
	spxy = new double[Nconfs];
	spyz = new double[Nconfs];
	
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
	
	Ipxy = new double[N];
	Ipxz = new double[N];
	Ipyz = new double[N];
	
	Idxy = new double[N];
	Idxz = new double[N];
	Idyz = new double[N];
	
	Irxy = new double[N];
	Irxz = new double[N];
	Iryz = new double[N];
	
	Time = new double[Nconfs];
		
	msd_pxy = new double*[N];
	msd_pxz = new double*[N];
	msd_pyz = new double*[N];
											
	msd_dxy = new double*[N];
	msd_dxz = new double*[N];
	msd_dyz = new double*[N];
											
	msd_rxy = new double*[N];
	msd_rxz = new double*[N];
	msd_ryz = new double*[N];
	
	msdp = new double*[N];
	msdd = new double*[N];
	msdr = new double*[N];
	
	msd = new double*[N];

// initialization of arrays 
	for(int i = 0; i < Nconfs; i++){
		
		if(i<N){
		msd_pxy[i] = new double[2];
		msd_pxz[i] = new double[2];
		msd_pyz[i] = new double[2];
		
		msd_dxy[i] = new double[2];
		msd_dxz[i] = new double[2];
		msd_dyz[i] = new double[2];
	
		msd_rxy[i] = new double[2];
		msd_rxz[i] = new double[2];
		msd_ryz[i] = new double[2];

		msdp[i] = new double[2];
		msdd[i] = new double[2];
		msdr[i] = new double[2];
		msd[i] = new double[2];
		
		Ipxy[i] = 0.0;
		Ipxz[i] = 0.0;
		Ipyz[i] = 0.0;
		
		Idxy[i] = 0.0;
		Idxz[i] = 0.0;
		Idyz[i] = 0.0;
		
		Irxy[i] = 0.0;
		Irxz[i] = 0.0;
		Iryz[i] = 0.0;
	}
		
		spxy[i] = 0.0;
		spxz[i] = 0.0;
		spyz[i] = 0.0;
		
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
		// read stress files
		ReadStressFile(input_file1,scxy,scxz,scyz,Time,Nconfs);
		ReadStressFile(input_file2,sdxy,sdxz,sdyz,Time,Nconfs);
		ReadStressFile(input_file3,srxy,srxz,sryz,Time,Nconfs);
		ReadStressFile(input_file4,skxy,skxz,skyz,Time,Nconfs);
		
		double dt = Time[1] - Time[0];
		
		cout << "I have read all stress files!" << endl;
		for(int j = 0; j<Nconfs; j++){
			spxy[j] = scxy[j] + skxy[j] ;
			spxz[j] = scxz[j] + skxz[j] ;
			spyz[j] = scyz[j] + skyz[j] ;
		}
		   
	// find integrals of stress tensor. 3 integrals are output for potential, dissipative and random stresses
	trapz(Ipxy,Ipxz, Ipyz, spxy,spxz,spyz,Idxy,Idxz,Idyz,sdxy,sdxz,sdyz,Irxy,Irxz,Iryz,srxy,srxz,sryz,Time,N);

	cout << "found all integrals of stresses " << endl;
	
	//The diffusion of momentum is calculated for the 3 directions and the 3 force contributions.
	//The square of the integral is calculated first and then the displacement is found.	
	momentum_msd(Ipxy,Ipxy,msd_pxy, N);
	momentum_msd(Ipxz,Ipxz,msd_pxz, N);
	momentum_msd(Ipyz,Ipyz,msd_pyz, N);
	
	cout << "found momentum diffusion for potential force" << endl;
	
	momentum_msd(Idxy,Idxy,msd_dxy, N);
	momentum_msd(Idxz,Idxz,msd_dxz, N);
	momentum_msd(Idyz,Idyz,msd_dyz, N);
	
	cout << "found momentum diffusion for dissipative force" << endl;
	
	momentum_msd(Irxy,Irxy,msd_rxy, N);
	momentum_msd(Irxz,Irxz,msd_rxz, N);
	momentum_msd(Iryz,Iryz,msd_ryz, N);
	
	cout << "found momentum diffusion for random force" << endl;

	double sum_r = 0.0,step;
	double slope1 = 0.0;

	for(int k = 0; k<N-2;k++){
		//average in the 3 directions 
		msdp[k][0] = (msd_pxy[k][0] + msd_pxz[k][0] + msd_pyz[k][0])/3.0;
		msdp[k][1] = msd_pxy[k][1];
		
		msdd[k][0] = (msd_dxy[k][0] + msd_dxz[k][0] + msd_dyz[k][0])/3.0;
		msdd[k][1] = msd_dxy[k][1];
		
		msdr[k][0] = (msd_rxy[k][0] + msd_rxz[k][0] + msd_ryz[k][0])/3.0;
		msdr[k][1] = msd_rxy[k][1];
		
		// eta_pp - eta_dd is calculated according to the reference M. Panoukidou, CR Wand, P. Carbobe, Phys.Rew E (2020)
		msd[k][0] = msdp[k][0] - msdd[k][0];
		msd[k][1] = msdp[k][1];

			
		// find the slope of the momentum diffusion using finite differences 
		if(k>0 && k<N-1){
			slope1 = (msd[k][0] - msd[k-1][0])/(Time[k] - Time[k-1]);
		}else{
			slope1 = 0.0;
		}	
		slope1 = slope1*V/(2.0) + msdr[0][0]*V/dt;
		
		writeMSDFile(msd_file,msd,slope1,msdr[k][0],k);
	}

	//dealocate memory
	for (int i = 0; i < N; i++){

		delete [] msd_pxy [i];
		delete [] msd_pxz [i];
		delete [] msd_pyz [i];

		delete [] msd_dxy [i];
		delete [] msd_dxz [i];
		delete [] msd_dyz [i];

		delete [] msd_rxy [i];
		delete [] msd_rxz [i];
		delete [] msd_ryz [i];

		delete [] msdp [i];
		delete [] msdd [i];	
		delete [] msdr [i];	
		
		delete [] msd [i];

	}
	delete [] msd_pxy;
	delete [] msd_pxz;
	delete [] msd_pyz;
	delete [] msd_dxy;
	delete [] msd_dxz;
	delete [] msd_dyz;
	delete [] msd_rxy;
	delete [] msd_rxz;
	delete [] msd_ryz;
	delete [] msd;
	delete [] msdp;
	delete [] msdd;
	delete [] msdr;
	delete [] Time;
	delete [] Ipxy;
	delete [] Ipxz;
	delete [] Ipyz;	
	delete [] Idxy;
	delete [] Idxz;
	delete [] Idyz;	
	delete [] Irxy;
	delete [] Irxz;
	delete [] Iryz;
	delete [] spxz;
	delete [] spyz;
	delete [] scxy;
	delete [] spxy;
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
