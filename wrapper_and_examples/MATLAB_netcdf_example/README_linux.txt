If have trouble with Octave/MATLAB reading netcdf files ensure that the dynamic link loader path points to the correct libnetcdf.so file
e.g.
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
