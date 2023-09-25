addpath("..")
pkg load netcdf

% read in data from a netcdf file
target=ncread("test5_input_data.nc","target");
timebase=ncread("test5_input_data.nc","timebase");
p1=ncread("test5_input_data.nc","proxy1");
p2=ncread("test5_input_data.nc","proxy2");

% call SLICKER with default settings except for number of ensemble members, time limit for each ensemble member and the solution tolerance
[time,central_tendency,uncertainty,spread]=SLICKER(target,timebase,p1,p2,'num_ensemble',512,'time_limit',20,'tol',1e-7);

% write the output to file
t=[time,central_tendency,uncertainty,spread]';
fileID = fopen('test5_MATLAB_output.txt','w');
fprintf(fileID,'%12s %12s %12s %8s\n','!time base','M-estimator','+/- for_95%_CI',' Qn');
fprintf(fileID,'%12.4e %12.4e %12.4e %12.4e\n',t);
fclose(fileID);


