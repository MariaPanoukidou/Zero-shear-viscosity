#ifndef EIN_LIB
#define EIN_LIB

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "ein_lib.h"

using namespace std;
extern string in_stress, visc,msd_file;
extern ifstream infile;
extern ofstream MSDF;
extern double *Time;
extern double *scxy, *scxz, *scyz;
extern double *sdxy, *sdxz, *sdyz;
extern double *srxy, *srxz, *sryz;
extern double *skxy, *skxz, *skyz;
extern double *Ipxy,*Ipxz,*Ipyz,*Idxy,*Idxz,*Idyz,*Irxy ,*Irxz, *Iryz;
extern double **msd,**msdp,**msdd, **msdr;
extern double **msd_pxy,**msd_pxz, **msd_pyz, **msd_dxy, **msd_dxz,**msd_dyz, **msd_rxy, **msd_rxz,**msd_ryz;
extern int Nconfs;

void ReadStressFile(string in_stress, double *sxy,double *sxz, double *syz, double *t, int Nconfs);
void writeMSDFile(string msd_file, double **msd, double slope1, double msdr, int t);
void momentum_msd(double *y,double *x,double **msd, int Nconfs);
void trapz(double *I1xy,double *I1xz,double *I1yz, double *sxy1,double *sxz1, double *syz1,double *I2xy,double *I2xz, double *I2yz, double *sxy2,double *sxz2, double *syz2,double *I3xy,double *I3xz,double *I3yz, double *sxy3,double *sxz3, double *syz3,double *Time, int N);

#endif

