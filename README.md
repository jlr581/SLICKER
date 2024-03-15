# Segmented Linear Integral Correlation Kernel Ensemble Reconstruction

## Project Description
Segmented Linear Integral Correlation Kernel Ensemble Reconstruction (SLICKER) is a numerical method to reconstruct a (typically climate) time-series from multiple proxy data-series.  SLICKER allows for
- Univariate or multivariate reconstructions 
- Uneven and different data-spacing between allow proxies and the reconstruction target
- Different linear and non-linear relationships between the target and each proxy
- Ensembles to allow for robust estimates of the central tendency, uncertainty in the central tendency, and spread of the reconstruction
- Parallel execution is available via OpenMP

SLICKER is written in FORTRAN, with wrappers for MATLAB, Python and R.  A windows executable is also supplied for users without access to a FORTRAN compiler.

## Citation
If using SLICKER for scientific studies, please cite

Roberts et al (submitted), New insights from an East Antarctic ice core with SLICKER: Segmented Linear Integral Correlation Kernel Ensemble Reconstruction, *PLOS ONE*

## Building Source Code
For systems with access to a FORTRAN compiler, source code is included
in the ../FORTRAN_source_code directory.  
 - Edit the "system variables" section of the Makefile to set system 
   specific options, noting that defaults for ifort, gfortran and 
   nvfortran are provided, for both serial and OpenMP parallel
   configurations.  
 - Compile the code by typing "make" at the command line, will generate the executable "reconstruct"
 - copy "reconstruct" to the "wrappers_and_examples" directory by typing "make install" at the command line

## Wrappers and examples
Wrappers and example usage is provided for MATLAB, Python and R.  The particular text example for the wrappers is from Figure 8c of the above *The Holocene* journal paper a test case with two input proxies with 20% of the data missing from both proxies. Note that the solver relies on random numbers, so the solution from subsequent runs of the same problem will be slightly different.  Examples for both text and netCDF input data are provided.  All other test examples from *The Holocene* journal paper are included in the Journal_paper_test_cases directory.

| Directory | code execution |
|-----------|----------------|
| MATLAB_netcdf_example | matlab -nodisplay test5.m |
| MATLAB_text_example | matlab -nodisplay test5.m |
| Python | python reconstruct.py |
| Python | jupyter lab reconstruct_example.ipynb |
| R | R --no-save < test5.R |

