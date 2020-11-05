# Zero-shear-viscosity
Green Kubo and Einstein-Helfand methods for calculation of viscosity.
The codes are developed to post-process results output from DL_MESO_DPD software. 

input files: Stess_pot.d
             Stress_kin.d
             Stress_rn.d
             Stress_diss.d
             
output files of GK code: pressure_autocorr.dat
                         viscosity_integral.dat
                         
output files of Einstein code: msd_mom.dat

To compile GK:
              g++ gk.cpp gk_lib.h gk_lib.cpp -o gk.exe
              
Similarly, to compile Einstein-Helfand method:
              g++ visc_einstein.cpp ein_lib.h ein_lib.cpp -o ein.exe
