#ifndef GK_LIB
#define GK_LIB

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
#include "gk_lib.h"

#define PI 3.14159265
using namespace std;
extern string in_stress, stress_d, stress_k, stress_r, AC, visc;
extern ifstream infile, infile_d, infile_r, inflie_k;
extern ofstream autocorrF, viscF;
extern double *Time;
extern double *scxy, *scxz, *scyz;
extern double *sdxy, *sdxz, *sdyz;
extern double *srxy, *srxz, *sryz;
extern double *skxy, *skxz, *skyz;
extern double AcR0_xy,AcR0_xz,AcR0_yz;
extern int Nconfs;

void ReadStressFile(string in_stress, double *sxy,double *sxz, double *syz, double *t, int Nconfs);
void writeAutocorrFile(string AC, double *mtoxy,double *mtoxz, double *mtoyz, double *tt, int t);
void writeViscFile(string visc, double etaxy,double etaxz,double etayz, double *tt, int t);
double CalcAutocorr(double *mtoValue, double *s, double *sd, double *sr, int i, int Nconfs);

#endif

