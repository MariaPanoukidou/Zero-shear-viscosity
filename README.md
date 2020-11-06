# Zero-shear-viscosity
Green Kubo and Einstein-Helfand methods for calculation of viscosity.
The codes are developed to post-process results output from DL_MESO_DPD software. 

=================== VERY IMPORTANT NOTE ===========================================

Note: Before running the code please make sure that you have deleted the output files from previous runs since the code will append the new data at the end of the old file(s) and this might create problems in the analysis of the results. 

==================== input and output files =======================================

input files: Stess_pot.d
             Stress_kin.d
             Stress_rn.d
             Stress_diss.d
        These file are output from DL_MESO(DPD). 
        
output files of GK code: pressure_autocorr.dat
                         viscosity_integral.dat
                         
output files of Einstein code: momentum_diff.dat

================= technical details ===============================================

To compile GK:
              g++ gk.cpp gk_lib.h gk_lib.cpp -o gk.exe
              
Similarly, to compile Einstein-Helfand method:
              g++ visc_einstein.cpp ein_lib.h ein_lib.cpp -o ein.exe

To run: ./gk.exe <Nconfs> <Volume> <kbT> <Nintegral>
                or
 ./ein.exe <Nconfs> <Volume> <kbT> <Nintegral>
  
where Nconfs = number of lines in the Stress_* files
      Volume = the volume of the simulation box
      kbT = the target thermal energy used for the DPD simulations
      Nintegral = the number of points the user would like to include in the integral calculations. (Nintegral <= Nconfs)

==================== information on the codes======================================

gk.cpp: The main function where the generalized Green - Kubo formula is applied.
gk_lib.h: header file with function and variable declarations.
gk_lib.cpp: functions called in the main code (gk.cpp). These functions are the reading and writing of files, calculation of stress autocorrelation function and integral using trapezoidal rule. 

visc_einstein.cpp: The main function where the generalized Einstein - Helfand formula is applied.
ein_lib.h: header file with function and variable declarations.
ein_lib.cpp: functions called in the main code (visc_einstein.cpp). These functions are the reading and writing of files, calculation of stress integrals using trapezoidal rule and calculation of momentum diffusion.
